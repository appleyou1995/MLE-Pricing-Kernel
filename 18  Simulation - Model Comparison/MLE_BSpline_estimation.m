%% Main Function: MLE Theta Estimation for B-Spline Pricing Kernel

function [theta_hat, log_lik, BIC, exitflag, output, delta_vec, M_vec, pit_vec] = MLE_BSpline_estimation( ...
    Smooth_AllR, Smooth_AllR_RND, Realized_Return, Risk_Free_Rate, ...
    b, alpha, beta, Global_Min_R, Global_Max_R, ...
    impose_monotonicity, Basis_Precomputed_Input)

    if nargin < 10 || isempty(impose_monotonicity)
        impose_monotonicity = true;
    end

    if nargin < 11
        Basis_Precomputed_Input = [];
    end

    % If beta contains more than one value, treat it as a grid
    % and select the beta with the largest maximized likelihood.
    beta_candidates = beta(:);

    if isempty(beta_candidates) || ...
            any(~isfinite(beta_candidates)) || ...
            any(beta_candidates <= 0)

        error('beta must contain finite positive values.');
    end

    if numel(beta_candidates) > 1
        num_beta = numel(beta_candidates);

        theta_grid = cell(num_beta, 1);
        output_grid = cell(num_beta, 1);

        ll_grid = -Inf(num_beta, 1);
        bic_grid = Inf(num_beta, 1);
        exit_grid = -999 * ones(num_beta, 1);

        for beta_idx = 1:num_beta
            beta_now = beta_candidates(beta_idx);

            [theta_now, ll_now, bic_now, ...
                exit_now, output_now] = ...
                MLE_BSpline_estimation( ...
                Smooth_AllR, Smooth_AllR_RND, ...
                Realized_Return, Risk_Free_Rate, ...
                b, alpha, beta_now, ...
                Global_Min_R, Global_Max_R, ...
                impose_monotonicity, ...
                Basis_Precomputed_Input);

            theta_grid{beta_idx} = theta_now;
            output_grid{beta_idx} = output_now;

            ll_grid(beta_idx) = ll_now;
            bic_grid(beta_idx) = bic_now;
            exit_grid(beta_idx) = exit_now;
        end

        valid = isfinite(ll_grid) & exit_grid > 0;

        if ~any(valid)
            theta_hat = zeros(b + 1, 1);
            log_lik = -Inf;
            BIC = Inf;
            exitflag = -999;

            output = struct( ...
                'SelectedBeta', NaN, ...
                'BetaGrid', beta_candidates, ...
                'LogLikelihoodGrid', ll_grid, ...
                'BICGrid', bic_grid, ...
                'ExitFlagGrid', exit_grid);

            delta_vec = [];
            M_vec = [];
            pit_vec = [];
            return
        end

        ll_for_selection = ll_grid;
        ll_for_selection(~valid) = -Inf;

        [log_lik, best_idx] = max(ll_for_selection);

        theta_hat = theta_grid{best_idx};
        exitflag = exit_grid(best_idx);
        beta_hat = beta_candidates(best_idx);
        output = output_grid{best_idx};

        output.SelectedBeta = beta_hat;
        output.BetaGrid = beta_candidates;
        output.LogLikelihoodGrid = ll_grid;
        output.BICGrid = bic_grid;
        output.ExitFlagGrid = exit_grid;

        % beta is estimated in this branch, so count it in BIC.
        BIC = (length(theta_hat) + 1) * ...
            log(length(Realized_Return.realized_ret)) ...
            - 2 * log_lik;

        if nargout > 5
            [~, ~, delta_vec, M_vec, pit_vec] = ...
                log_likelihood_bspline( ...
                theta_hat, ...
                Realized_Return.realized_ret(:), ...
                Risk_Free_Rate(:), ...
                Basis_Precomputed_Input, ...
                Smooth_AllR, ...
                Smooth_AllR_RND, ...
                Smooth_AllR.Properties.VariableNames, ...
                alpha, beta_hat);
        else
            delta_vec = [];
            M_vec = [];
            pit_vec = [];
        end

        return
    end

    % Preserve the original behavior when beta is scalar.
    beta = beta_candidates(1);

    % Settings
    rng(0);
    R_vec  = Realized_Return.realized_ret;
    Rf_vec = Risk_Free_Rate;
    months = Smooth_AllR.Properties.VariableNames;
    T      = length(R_vec);

    % --- Step 1: Construct knots and precompute basis ---
    n_degree = 5;
    k_order  = n_degree + 1;
    min_knot = Global_Min_R;
    max_knot = Global_Max_R;

    num_basis_function = b + 1;
    num_breaks = num_basis_function - k_order + 2;

    if num_breaks < 2
        error(['The number of basis functions is insufficient: ', ...
            'b=%d, degree=%d.'], b, n_degree);
    end

    breaks = linspace(min_knot, max_knot, num_breaks);
    knots  = augknt(breaks, k_order);

    if isempty(Basis_Precomputed_Input)
        Basis_Precomputed = cell(T, 1);

        for t = 1:T
            try
                col_name = months{t};
                R_axis_t = Smooth_AllR.(col_name);
                Basis_Precomputed{t} = spcol( ...
                    knots, k_order, R_axis_t(:));
            catch ME
                warning('Error in basis construction at t=%d (%s): %s', ...
                    t, months{t}, ME.message);
                Basis_Precomputed{t} = [];
            end
        end
    else
        Basis_Precomputed = Basis_Precomputed_Input;

        if numel(Basis_Precomputed) ~= T
            error(['Basis_Precomputed_Input must contain one cell for ', ...
                'each of the %d months.'], T);
        end
    end

    % --- Step 2: Optional monotonicity restriction ---
    num_params = b + 1;

    if impose_monotonicity
        % theta_{i+1} - theta_i <= 0
        A_ineq = zeros(b, num_params);
        for i = 1:b
            A_ineq(i, i)   = -1;
            A_ineq(i, i+1) =  1;
        end
        b_ineq = zeros(b, 1);
    else
        A_ineq = [];
        b_ineq = [];
    end

    % --- Step 3: Optimization setup ---
    theta0 = zeros(num_params, 1);

    options = optimoptions('fmincon', ...
        'Display', 'off', ...
        'Algorithm', 'sqp', ...
        'ConstraintTolerance', 1e-6, ...
        'StepTolerance', 1e-6, ...
        'MaxFunctionEvaluations', 50000);

    obj_fun = @(param) -log_likelihood_bspline( ...
        param, R_vec, Rf_vec, ...
        Basis_Precomputed, Smooth_AllR, Smooth_AllR_RND, ...
        months, alpha, beta);

    % --- Step 4: Run optimization ---
    try
        [theta_hat, neg_LL, exitflag, output] = fmincon( ...
            obj_fun, theta0, ...
            A_ineq, b_ineq, [], [], [], [], [], options);
    catch ME
        warning(['fmincon failed for b=%d, alpha=%.2f, beta=%.2f, ', ...
            'monotonicity=%d. Message: %s'], ...
            b, alpha, beta, impose_monotonicity, ME.message);

        theta_hat = theta0;
        log_lik   = -Inf;
        BIC       = Inf;
        exitflag  = -999;
        output    = struct();
        output.SelectedBeta = beta;

        delta_vec = [];
        M_vec     = [];
        pit_vec   = [];
        return
    end

    log_lik = -neg_LL;
    output.SelectedBeta = beta;

    % Post-estimation call
    if nargout > 5
        [~, BIC, delta_vec, M_vec, pit_vec] = log_likelihood_bspline( ...
            theta_hat, R_vec, Rf_vec, ...
            Basis_Precomputed, Smooth_AllR, Smooth_AllR_RND, ...
            months, alpha, beta);
    else
        [~, BIC] = log_likelihood_bspline( ...
            theta_hat, R_vec, Rf_vec, ...
            Basis_Precomputed, Smooth_AllR, Smooth_AllR_RND, ...
            months, alpha, beta);

        delta_vec = [];
        M_vec     = [];
        pit_vec   = [];
    end
end
