%% Main Function: GMM Theta Estimation for Pricing Kernel

function [theta_hat, Jval, delta_vec, num_vec, den_vec] = GMM_theta_estimation( ...
    Smooth_AllR, Smooth_AllR_RND, Realized_Return, Risk_Free_Rate, b, min_knot, max_knot)

    % Set the seed for the random number generator
    rng(0);

    % Initial parameters
    theta0 = zeros(1, b+1);

    LB = [];
    UB = [];

    % Pull inputs
    dates   = Realized_Return.date;
    R_vec  = Realized_Return.realized_ret;                                 % realized gross return R_{t+1}
    Rf_vec = Risk_Free_Rate.rate;                                          % risk-free R^f_t

    % Optimization options
    options = optimoptions('fmincon', ...
        'Display', 'iter-detailed', ...
        'Algorithm', 'interior-point', ...
        'SpecifyObjectiveGradient', false);

    % Define objective function
    obj_fun = @(theta) GMM_objective_function(theta, ...
                     Smooth_AllR, Smooth_AllR_RND, R_vec, Rf_vec, dates, b, min_knot, max_knot);

    % Run optimization
    [theta_hat, Jval, exitflag] = fmincon(obj_fun, theta0, ...
                                [],[],[],[], LB, UB, [], options);

    % Post-estimation call
    [~, delta_vec, num_vec, den_vec] = GMM_moment_conditions(theta_hat, ...
        Smooth_AllR, Smooth_AllR_RND, R_vec, Rf_vec, dates, b, min_knot, max_knot);

    % Report
    if     exitflag > 0,  disp('✅ fmincon converged successfully.');
    elseif exitflag == 0, disp('⚠️ fmincon reached the iteration limit.');
    else,                 disp(['❌ fmincon failed. Exit flag: ', num2str(exitflag)]);
    end
    
    % Display output
    disp('Estimated theta:'); disp(theta_hat);
    disp(['b = ', num2str(b)]);

end


%% Local Function: GMM Objective Function

function J = GMM_objective_function(theta, ...
    Smooth_AllR, Smooth_AllR_RND, R_vec, Rf_vec, dates, b, min_knot, max_knot)

    [g, ~, ~, ~] = GMM_moment_conditions(theta, ...
        Smooth_AllR, Smooth_AllR_RND, R_vec, Rf_vec, dates, b, min_knot, max_knot);

    % Use a GMM type optimization with only the first stage optimization
    W = eye(b + 1);

    % Objective function
    J = g' * W * g;
end


%% Local Function: GMM Moment Conditions

function [g, delta_vec, num_vec, den_vec] = GMM_moment_conditions(theta, ...
    Smooth_AllR, Smooth_AllR_RND, R_vec, Rf_vec, dates, b, min_knot, max_knot)
    
    T = numel(R_vec);
    n_degree = 3;                                                          % Cubic B-spline
    m = b;
    g = zeros(m + 1, 1);

    delta_vec = NaN(T, 1);
    num_vec   = NaN(T, 1);
    den_vec   = NaN(T, 1);

    % Accumulators for sample moments
    moment_sum = zeros(m + 1, 1);
    valid_T    = zeros(m + 1, 1);

    for t = 1:T
        try
            R_realized_t = R_vec(t);
            Rf_t = Rf_vec(t);

            % Extract monthly return grid and RND values
            date_str = num2str(dates(t));
            R_grid = Smooth_AllR.(date_str);
            fstar  = Smooth_AllR_RND.(date_str);

            % Build B-spline basis matrix (L x (b+1))
            N = numel(R_grid);
            Bmat = zeros(b+1, N);

            for i = 1:(b + 1)
                Bmat(i, :) = Bspline_basis_function_value(n_degree, b, min_knot, max_knot, i, R_grid);
            end

            % phi(R) = sum_i theta_i * B_i(R)
            phi = theta * Bmat;

            % Log-sum-exp stabilization
            z = -phi;
            c = max(z);
            w = exp(z - c) .* fstar;                                       % weights proportional to exp(-phi) * f*_t(R)
            I_minus = trapz(R_grid, w) * exp(c);
            if ~(I_minus > 0 && isfinite(I_minus)), continue; end

            % δ_t = -ln Rf_t + ln I_minus
            delta_vec(t) = -log(Rf_t) + log(I_minus);

            % Pricing kernel on the grid:  M(R) = exp{ delta_t + phi(R) }
            M_grid = exp(delta_vec(t) + phi);

            % integrand: (R_f * M)^(-1) * f^Q
            integrand = fstar ./ (Rf_t * M_grid);

            % U ∈ [0,1]
            idx = (R_grid <= R_realized_t);
            num = trapz(R_grid(idx), integrand(idx));
            den = trapz(R_grid,      integrand);
            U = num / den;
            U = max(0, min(1, U));

            num_vec(t) = num;
            den_vec(t) = den;

            % Accumulate moments
            for j = 0:m
                moment_sum(j+1) = moment_sum(j+1) + U^(j+1);
                valid_T(j+1)    = valid_T(j+1)    + 1;
            end

        catch
            % Skip this month on error
            continue;
        end
    end

    % Sample moment vector g(theta)
    for j = 0:m
        if valid_T(j+1) > 0
            g(j+1) = (moment_sum(j+1) / valid_T(j+1)) - 1 / (j + 2);
        else
            % Penalize missing moments to avoid NaN
            g(j+1) = 1e3;
        end
    end
end
