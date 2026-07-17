%% Main Function: MLE Theta Estimation for B-Spline Pricing Kernel
%
% This version is backward-compatible with the original function.
%
%   Basis_Precomputed_Input
%       Optional precomputed B-spline basis. Supplying it avoids rebuilding
%       the same basis in every simulation replication.

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

        delta_vec = [];
        M_vec     = [];
        pit_vec   = [];
        return
    end

    log_lik = -neg_LL;

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
