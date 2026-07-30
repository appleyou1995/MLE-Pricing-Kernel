%% Main Function: MLE Theta Estimation for B-Spline Pricing Kernel

function [theta_hat, log_lik, BIC, exitflag, output, ...
    delta_vec, M_vec, pit_vec] = MLE_BSpline_estimation( ...
    Smooth_AllR, Smooth_AllR_RND, Realized_Return, Risk_Free_Rate, ...
    b, alpha, beta, Global_Min_R, Global_Max_R, ...
    impose_monotonicity, Basis_Precomputed_Input, ...
    Search_Settings, Monthly_Precomputed_Input)

    if nargin < 10 || isempty(impose_monotonicity)
        impose_monotonicity = true;
    end

    if nargin < 11
        Basis_Precomputed_Input = [];
    end

    if nargin < 12
        Search_Settings = [];
    end

    if nargin < 13
        Monthly_Precomputed_Input = [];
    end

    if ~isempty(Search_Settings)
        [theta_hat, log_lik, BIC, exitflag, output, ...
            delta_vec, M_vec, pit_vec] = ...
            estimate_with_parameter_search( ...
                Smooth_AllR, Smooth_AllR_RND, ...
                Realized_Return, Risk_Free_Rate, ...
                b, alpha, beta, ...
                Global_Min_R, Global_Max_R, ...
                impose_monotonicity, ...
                Basis_Precomputed_Input, ...
                Search_Settings, ...
                Monthly_Precomputed_Input, nargout);
        return
    end

    if ~isscalar(alpha) || ~isfinite(alpha) || alpha <= 0
        error('alpha must be a finite positive scalar.');
    end

    if ~isscalar(beta) || ~isfinite(beta) || beta <= 0
        error('beta must be a finite positive scalar.');
    end

    R_vec  = Realized_Return.realized_ret(:);
    Rf_vec = Risk_Free_Rate(:);
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
        months, alpha, beta, Monthly_Precomputed_Input);

    % --- Step 4: Run optimization ---
    try
        [theta_hat, neg_LL, exitflag, output] = fmincon( ...
            obj_fun, theta0, ...
            A_ineq, b_ineq, [], [], [], [], [], options);
    catch ME
        warning(['fmincon failed for b=%d, alpha=%.3f, beta=%.3f, ', ...
            'monotonicity=%d. Message: %s'], ...
            b, alpha, beta, impose_monotonicity, ME.message);

        theta_hat = theta0;
        log_lik   = -Inf;
        BIC       = Inf;
        exitflag  = -999;
        output    = struct();
        output.SelectedAlpha = alpha;
        output.SelectedBeta = beta;
        output.SearchParameter = "none";

        delta_vec = [];
        M_vec     = [];
        pit_vec   = [];
        return
    end

    log_lik = -neg_LL;
    output.SelectedAlpha = alpha;
    output.SelectedBeta = beta;
    output.SearchParameter = "none";
    output.ParameterGrid = [];
    output.LogLikelihoodGrid = [];
    output.ExitFlagGrid = [];
    output.SearchExpanded = false;
    output.BoundaryHit = false;
    output.SearchStopReason = "fixed parameter";

    % Post-estimation call
    if nargout > 5
        [~, BIC, delta_vec, M_vec, pit_vec] = ...
            log_likelihood_bspline( ...
                theta_hat, R_vec, Rf_vec, ...
                Basis_Precomputed, ...
                Smooth_AllR, Smooth_AllR_RND, ...
                months, alpha, beta, ...
                Monthly_Precomputed_Input);
    else
        [~, BIC] = log_likelihood_bspline( ...
            theta_hat, R_vec, Rf_vec, ...
            Basis_Precomputed, ...
            Smooth_AllR, Smooth_AllR_RND, ...
            months, alpha, beta, ...
            Monthly_Precomputed_Input);

        delta_vec = [];
        M_vec     = [];
        pit_vec   = [];
    end
end


function [theta_hat, log_lik, BIC, exitflag, output, ...
    delta_vec, M_vec, pit_vec] = estimate_with_parameter_search( ...
    Smooth_AllR, Smooth_AllR_RND, ...
    Realized_Return, Risk_Free_Rate, ...
    b, alpha_fixed, beta_fixed, ...
    Global_Min_R, Global_Max_R, ...
    impose_monotonicity, Basis_Precomputed_Input, ...
    S, Monthly_Precomputed_Input, caller_nargout)

    required_fields = { ...
        'Parameter', 'NullValue', 'CoarseGrid', ...
        'CoarseStep', 'FineStep', ...
        'LowerBoundExclusive', 'UpperBoundInclusive', ...
        'LikelihoodTolerance', ...
        'NonImprovingToStop', 'MaxExpansionSteps'};

    if ~all(isfield(S, required_fields))
        error('Search_Settings is missing required fields.');
    end

    parameter = lower(string(S.Parameter));
    if parameter ~= "alpha" && parameter ~= "beta"
        error('Search parameter must be alpha or beta.');
    end

    lower_bound = S.LowerBoundExclusive;
    upper_bound = S.UpperBoundInclusive;
    tolerance = S.LikelihoodTolerance;

    if ~isscalar(lower_bound) || ~isscalar(upper_bound) || ...
            lower_bound < 0 || upper_bound <= lower_bound
        error('Invalid parameter-search bounds.');
    end

    candidates = unique(round(S.CoarseGrid(:), 10), 'sorted');
    candidates = candidates( ...
        candidates > lower_bound & candidates <= upper_bound);

    if isempty(candidates)
        error('No coarse-grid candidates lie within the search bounds.');
    end

    null_value = S.NullValue;
    if null_value <= lower_bound || null_value > upper_bound
        error('The null value lies outside the search bounds.');
    end

    if ~any(abs(candidates - null_value) <= 1e-10)
        candidates = unique([candidates; null_value], 'sorted');
    end

    Show_Search_Progress = isempty(getCurrentTask());
    Search_Timer = tic;

    evaluated_values = zeros(0, 1);
    theta_grid = cell(0, 1);
    optimizer_output_grid = cell(0, 1);
    ll_grid = zeros(0, 1);
    bic_grid = zeros(0, 1);
    exit_grid = zeros(0, 1);
    cached_candidate_used = false;

    % A fitted null model can be supplied as the alpha=1 (or beta=1)
    % candidate. This avoids solving an identical constrained MLE again
    % inside every parameter search.
    if isfield(S, 'CachedCandidate') && ...
            ~isempty(S.CachedCandidate)

        C = S.CachedCandidate;
        cached_fields = { ...
            'Parameter', 'Value', 'Theta', ...
            'LogLikelihood', 'BIC', 'ExitFlag', ...
            'OptimizerOutput'};

        if ~all(isfield(C, cached_fields))
            error('CachedCandidate is missing required fields.');
        end

        cached_parameter = lower(string(C.Parameter));
        cached_value = C.Value;

        if cached_parameter ~= parameter
            error('CachedCandidate parameter does not match the search.');
        end
        if ~isscalar(cached_value) || ~isfinite(cached_value) || ...
                cached_value <= lower_bound || ...
                cached_value > upper_bound
            error('CachedCandidate value lies outside the search bounds.');
        end
        if numel(C.Theta) ~= b + 1
            error('CachedCandidate theta has an unexpected length.');
        end
        if ~isscalar(C.LogLikelihood) || ...
                ~isfinite(C.LogLikelihood) || ...
                ~isscalar(C.ExitFlag) || ...
                ~isfinite(C.ExitFlag)
            error('CachedCandidate has invalid likelihood or exitflag.');
        end
        if ~isstruct(C.OptimizerOutput)
            error('CachedCandidate OptimizerOutput must be a struct.');
        end

        evaluated_values = round(cached_value, 10);
        theta_grid = {C.Theta(:)};
        optimizer_output_grid = {C.OptimizerOutput};
        ll_grid = C.LogLikelihood;
        bic_grid = C.BIC;
        exit_grid = C.ExitFlag;
        cached_candidate_used = true;
    end

    if Show_Search_Progress
        fprintf('\nStarting adaptive %s search.\n', parameter);
    end

    [evaluated_values, theta_grid, optimizer_output_grid, ...
        ll_grid, bic_grid, exit_grid] = evaluate_candidates( ...
            candidates, evaluated_values, ...
            theta_grid, optimizer_output_grid, ...
            ll_grid, bic_grid, exit_grid, ...
            parameter, alpha_fixed, beta_fixed, ...
            Smooth_AllR, Smooth_AllR_RND, ...
            Realized_Return, Risk_Free_Rate, ...
            b, Global_Min_R, Global_Max_R, ...
            impose_monotonicity, Basis_Precomputed_Input, ...
            Monthly_Precomputed_Input, Show_Search_Progress, ...
            Search_Timer);

    initial_min = min(evaluated_values);
    initial_max = max(evaluated_values);
    boundary_expanded = false;
    boundary_hit = false;
    stop_reasons = strings(0, 1);

    % Expand any competitive coarse-grid boundary. Evaluated candidates
    % are cached, so expansion never re-runs an existing point.
    for direction = [-1, 1]
        [~, best_ll] = select_best_candidate( ...
            evaluated_values, ll_grid, exit_grid, ...
            tolerance, null_value);

        if direction < 0
            edge_value = min(evaluated_values);
        else
            edge_value = max(evaluated_values);
        end

        edge_idx = find(abs(evaluated_values - edge_value) < 1e-10, ...
            1, 'first');

        if isempty(edge_idx) || ...
                ~isfinite(ll_grid(edge_idx)) || ...
                ll_grid(edge_idx) < best_ll - tolerance
            continue
        end

        consecutive_non_improving = 0;
        expansion_steps = 0;

        while consecutive_non_improving < S.NonImprovingToStop
            expansion_steps = expansion_steps + 1;
            if expansion_steps > S.MaxExpansionSteps
                boundary_hit = true;
                stop_reasons(end + 1, 1) = ...
                    "coarse expansion safety limit"; %#ok<AGROW>
                break
            end

            new_value = round( ...
                edge_value + direction * S.CoarseStep, 10);

            if new_value <= lower_bound || new_value > upper_bound
                boundary_hit = true;
                stop_reasons(end + 1, 1) = ...
                    "coarse parameter bound"; %#ok<AGROW>
                break
            end

            best_ll_before = best_ll;
            [evaluated_values, theta_grid, ...
                optimizer_output_grid, ll_grid, bic_grid, exit_grid] = ...
                evaluate_candidates( ...
                    new_value, evaluated_values, ...
                    theta_grid, optimizer_output_grid, ...
                    ll_grid, bic_grid, exit_grid, ...
                    parameter, alpha_fixed, beta_fixed, ...
                    Smooth_AllR, Smooth_AllR_RND, ...
                    Realized_Return, Risk_Free_Rate, ...
                    b, Global_Min_R, Global_Max_R, ...
                    impose_monotonicity, ...
                    Basis_Precomputed_Input, ...
                    Monthly_Precomputed_Input, ...
                    Show_Search_Progress, Search_Timer);

            boundary_expanded = true;
            new_idx = find( ...
                abs(evaluated_values - new_value) < 1e-10, ...
                1, 'first');
            new_ll = ll_grid(new_idx);

            [~, best_ll] = select_best_candidate( ...
                evaluated_values, ll_grid, exit_grid, ...
                tolerance, null_value);

            if isfinite(new_ll) && ...
                    new_ll > best_ll_before + tolerance
                consecutive_non_improving = 0;
            else
                consecutive_non_improving = ...
                    consecutive_non_improving + 1;
            end

            edge_value = new_value;
        end
    end

    [coarse_best_idx, ~] = select_best_candidate( ...
        evaluated_values, ll_grid, exit_grid, ...
        tolerance, null_value);
    coarse_best_value = evaluated_values(coarse_best_idx);

    fine_candidates = ( ...
        coarse_best_value - S.CoarseStep : ...
        S.FineStep : ...
        coarse_best_value + S.CoarseStep)';
    fine_candidates = unique(round(fine_candidates, 10), 'sorted');
    fine_candidates = fine_candidates( ...
        fine_candidates > lower_bound & ...
        fine_candidates <= upper_bound);

    [evaluated_values, theta_grid, optimizer_output_grid, ...
        ll_grid, bic_grid, exit_grid] = evaluate_candidates( ...
            fine_candidates, evaluated_values, ...
            theta_grid, optimizer_output_grid, ...
            ll_grid, bic_grid, exit_grid, ...
            parameter, alpha_fixed, beta_fixed, ...
            Smooth_AllR, Smooth_AllR_RND, ...
            Realized_Return, Risk_Free_Rate, ...
            b, Global_Min_R, Global_Max_R, ...
            impose_monotonicity, Basis_Precomputed_Input, ...
            Monthly_Precomputed_Input, Show_Search_Progress, ...
            Search_Timer);

    % If the final maximum is still at the complete explored boundary,
    % continue using the fine step and the same stopping rule.
    for direction = [-1, 1]
        [~, best_ll] = select_best_candidate( ...
            evaluated_values, ll_grid, exit_grid, ...
            tolerance, null_value);

        if direction < 0
            edge_value = min(evaluated_values);
        else
            edge_value = max(evaluated_values);
        end

        edge_idx = find(abs(evaluated_values - edge_value) < 1e-10, ...
            1, 'first');

        if isempty(edge_idx) || ...
                ~isfinite(ll_grid(edge_idx)) || ...
                ll_grid(edge_idx) < best_ll - tolerance
            continue
        end

        consecutive_non_improving = 0;
        expansion_steps = 0;

        while consecutive_non_improving < S.NonImprovingToStop
            expansion_steps = expansion_steps + 1;
            if expansion_steps > S.MaxExpansionSteps
                boundary_hit = true;
                stop_reasons(end + 1, 1) = ...
                    "fine expansion safety limit"; %#ok<AGROW>
                break
            end

            new_value = round( ...
                edge_value + direction * S.FineStep, 10);

            if new_value <= lower_bound || new_value > upper_bound
                boundary_hit = true;
                stop_reasons(end + 1, 1) = ...
                    "fine parameter bound"; %#ok<AGROW>
                break
            end

            best_ll_before = best_ll;
            [evaluated_values, theta_grid, ...
                optimizer_output_grid, ll_grid, bic_grid, exit_grid] = ...
                evaluate_candidates( ...
                    new_value, evaluated_values, ...
                    theta_grid, optimizer_output_grid, ...
                    ll_grid, bic_grid, exit_grid, ...
                    parameter, alpha_fixed, beta_fixed, ...
                    Smooth_AllR, Smooth_AllR_RND, ...
                    Realized_Return, Risk_Free_Rate, ...
                    b, Global_Min_R, Global_Max_R, ...
                    impose_monotonicity, ...
                    Basis_Precomputed_Input, ...
                    Monthly_Precomputed_Input, ...
                    Show_Search_Progress, Search_Timer);

            boundary_expanded = true;
            new_idx = find( ...
                abs(evaluated_values - new_value) < 1e-10, ...
                1, 'first');
            new_ll = ll_grid(new_idx);

            [~, best_ll] = select_best_candidate( ...
                evaluated_values, ll_grid, exit_grid, ...
                tolerance, null_value);

            if isfinite(new_ll) && ...
                    new_ll > best_ll_before + tolerance
                consecutive_non_improving = 0;
            else
                consecutive_non_improving = ...
                    consecutive_non_improving + 1;
            end

            edge_value = new_value;
        end
    end

    [best_idx, ~] = select_best_candidate( ...
        evaluated_values, ll_grid, exit_grid, ...
        tolerance, null_value);

    if isempty(best_idx)
        theta_hat = zeros(b + 1, 1);
        log_lik = -Inf;
        BIC = Inf;
        exitflag = -999;
        output = struct( ...
            'SelectedAlpha', NaN, ...
            'SelectedBeta', NaN, ...
            'SearchParameter', parameter, ...
            'ParameterGrid', evaluated_values, ...
            'LogLikelihoodGrid', ll_grid, ...
            'ExitFlagGrid', exit_grid, ...
            'SearchExpanded', boundary_expanded, ...
            'BoundaryHit', boundary_hit, ...
            'SearchStopReason', "no finite-likelihood candidate", ...
            'NumCandidatesEvaluated', numel(evaluated_values), ...
            'NumCandidatesEstimated', ...
                numel(evaluated_values) - cached_candidate_used, ...
            'CachedCandidateUsed', cached_candidate_used);
        delta_vec = [];
        M_vec = [];
        pit_vec = [];
        return
    end

    [sorted_values, sort_idx] = sort(evaluated_values);
    sorted_ll = ll_grid(sort_idx);
    sorted_exit = exit_grid(sort_idx);

    theta_hat = theta_grid{best_idx};
    log_lik = ll_grid(best_idx);
    exitflag = exit_grid(best_idx);
    output = optimizer_output_grid{best_idx};

    selected_value = evaluated_values(best_idx);
    selected_alpha = alpha_fixed;
    selected_beta = beta_fixed;
    if parameter == "alpha"
        selected_alpha = selected_value;
    else
        selected_beta = selected_value;
    end

    output.SelectedAlpha = selected_alpha;
    output.SelectedBeta = selected_beta;
    output.SearchParameter = parameter;
    output.ParameterGrid = sorted_values;
    output.LogLikelihoodGrid = sorted_ll;
    output.ExitFlagGrid = sorted_exit;
    output.SearchExpanded = boundary_expanded || ...
        min(sorted_values) < initial_min || ...
        max(sorted_values) > initial_max;
    output.BoundaryHit = boundary_hit;
    output.NumCandidatesEvaluated = numel(sorted_values);
    output.NumCandidatesEstimated = ...
        numel(sorted_values) - cached_candidate_used;
    output.CachedCandidateUsed = cached_candidate_used;
    if isempty(stop_reasons)
        output.SearchStopReason = "interior optimum";
    else
        output.SearchStopReason = strjoin(unique(stop_reasons), "; ");
    end

    % One estimated distortion parameter is added to the BIC count.
    BIC = (length(theta_hat) + 1) * ...
        log(length(Realized_Return.realized_ret)) - 2 * log_lik;

    if Show_Search_Progress
        fprintf([ ...
            'Adaptive %s search completed: selected %.3f from ', ...
            '%d evaluated candidates in %.2f minutes.\n\n'], ...
            parameter, selected_value, numel(sorted_values), ...
            toc(Search_Timer) / 60);
    end

    if caller_nargout > 5
        [~, ~, delta_vec, M_vec, pit_vec] = ...
            log_likelihood_bspline( ...
                theta_hat, ...
                Realized_Return.realized_ret(:), ...
                Risk_Free_Rate(:), ...
                Basis_Precomputed_Input, ...
                Smooth_AllR, Smooth_AllR_RND, ...
                Smooth_AllR.Properties.VariableNames, ...
                selected_alpha, selected_beta, ...
                Monthly_Precomputed_Input);
    else
        delta_vec = [];
        M_vec = [];
        pit_vec = [];
    end
end


function [evaluated_values, theta_grid, optimizer_output_grid, ...
    ll_grid, bic_grid, exit_grid] = evaluate_candidates( ...
    new_values, evaluated_values, theta_grid, ...
    optimizer_output_grid, ll_grid, bic_grid, exit_grid, ...
    parameter, alpha_fixed, beta_fixed, ...
    Smooth_AllR, Smooth_AllR_RND, ...
    Realized_Return, Risk_Free_Rate, ...
    b, Global_Min_R, Global_Max_R, ...
    impose_monotonicity, Basis_Precomputed_Input, ...
    Monthly_Precomputed_Input, Show_Search_Progress, Search_Timer)

    new_values = unique(round(new_values(:), 10), 'stable');

    for value_idx = 1:numel(new_values)
        value_now = new_values(value_idx);
        if any(abs(evaluated_values - value_now) < 1e-10)
            continue
        end

        alpha_now = alpha_fixed;
        beta_now = beta_fixed;
        if parameter == "alpha"
            alpha_now = value_now;
        else
            beta_now = value_now;
        end

        [theta_now, ll_now, bic_now, exit_now, output_now] = ...
            MLE_BSpline_estimation( ...
                Smooth_AllR, Smooth_AllR_RND, ...
                Realized_Return, Risk_Free_Rate, ...
                b, alpha_now, beta_now, ...
                Global_Min_R, Global_Max_R, ...
                impose_monotonicity, ...
                Basis_Precomputed_Input, [], ...
                Monthly_Precomputed_Input);

        evaluated_values(end + 1, 1) = value_now; %#ok<AGROW>
        theta_grid{end + 1, 1} = theta_now; %#ok<AGROW>
        optimizer_output_grid{end + 1, 1} = output_now; %#ok<AGROW>
        ll_grid(end + 1, 1) = ll_now; %#ok<AGROW>
        bic_grid(end + 1, 1) = bic_now; %#ok<AGROW>
        exit_grid(end + 1, 1) = exit_now; %#ok<AGROW>

        if Show_Search_Progress
            fprintf([ ...
                '  %s = %.3f, LL = %.6f, exitflag = %d, ', ...
                'elapsed = %.2f minutes\n'], ...
                parameter, value_now, ll_now, exit_now, ...
                toc(Search_Timer) / 60);
        end
    end
end


function [best_idx, best_ll] = select_best_candidate( ...
    values, ll_values, exit_values, tolerance, null_value)

    if numel(exit_values) ~= numel(ll_values)
        error('Exit-flag and likelihood grids must have equal lengths.');
    end

    % Model selection follows the maximized finite likelihood, as in the
    % original grid search. Exit flags remain saved as convergence
    % diagnostics, but an exitflag of zero must not make the adaptive
    % search believe that no candidate has yet been evaluated; otherwise
    % it can unnecessarily expand across the complete hard-bound range.
    valid = isfinite(ll_values);
    if ~any(valid)
        best_idx = [];
        best_ll = -Inf;
        return
    end

    best_ll = max(ll_values(valid));
    eligible = find(valid & ll_values >= best_ll - tolerance);

    distance_to_null = abs(values(eligible) - null_value);
    minimum_distance = min(distance_to_null);
    eligible = eligible( ...
        distance_to_null <= minimum_distance + 1e-12);

    if numel(eligible) > 1
        [~, local_idx] = min(values(eligible));
        best_idx = eligible(local_idx);
    else
        best_idx = eligible(1);
    end
end
