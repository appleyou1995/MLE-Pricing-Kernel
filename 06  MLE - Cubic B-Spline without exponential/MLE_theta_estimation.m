%% Main Function: MLE Theta Estimation for Pricing Kernel

function [theta_hat, log_lik, delta_vec, M_vec] = MLE_theta_estimation( ...
    Smooth_AllR, Smooth_AllR_RND, Realized_Return, Risk_Free_Rate, b)

    % Settings
    rng(0);
    theta0 = zeros(1, b+1);                                                % initial guess
    LB = [];
    UB = [];
    epsPos = 1e-6;

    % Pull inputs
    dates   = Realized_Return.date;
    R_vec   = Realized_Return.realized_ret;                                % realized gross return R_{t+1}
    Rf_vec  = Risk_Free_Rate;                                              % risk-free R^f_t

    % Optimization options
    options = optimoptions('fmincon', ...
        'Display', 'iter-detailed', ...
        'Algorithm', 'interior-point', ...
        'SpecifyObjectiveGradient', false);

    % Define objective function (maximize LL -> minimize -LL)
    obj_fun = @(theta) -log_likelihood_spline(theta, R_vec, Rf_vec, ...
                     Smooth_AllR, Smooth_AllR_RND, dates, b);

    % Run optimization
    nonlcon = @(theta) nlcon_pos(theta, Rf_vec, Smooth_AllR, Smooth_AllR_RND, dates, b, epsPos);
    [theta_hat, neg_LL, exitflag] = fmincon(obj_fun, theta0, ...
                                [],[],[],[], LB, UB, nonlcon, options);

    % Return log-likelihood (positive value)
    log_lik = -neg_LL;

    % Report
    if     exitflag > 0,  disp('✅ fmincon converged successfully.');
    elseif exitflag == 0, disp('⚠️ fmincon reached the iteration limit.');
    else,                 disp(['❌ fmincon failed. Exit flag: ', num2str(exitflag)]);
    end
    
    % Display output
    disp('Estimated theta:'); disp(theta_hat);
    disp(['b = ', num2str(b)]);
    disp(['Final log-likelihood = ', num2str(log_lik)]);

    % Post-estimation call
    [~, delta_vec, M_vec] = log_likelihood_spline(theta_hat, R_vec, Rf_vec, ...
    Smooth_AllR, Smooth_AllR_RND, dates, b);

end


%% Local Function: Log-Likelihood Function

function [LL, delta_vec, M_vec] = log_likelihood_spline(theta, R_vec, Rf_vec, ...
    Smooth_AllR, Smooth_AllR_RND, dates, b)
    
    T = length(R_vec);
    LL = 0;

    R_axis_1 = Smooth_AllR.(num2str(dates(1)));
    N        = numel(R_axis_1);

    delta_vec = NaN(T, 1);
    M_vec     = NaN(T, N);

    % log-likelihood
    for t = 1:T

        % === Step 1: Basic inputs ===
        R_realized_t = R_vec(t);                                           % (1*1) Realized gross return R_{t+1}
        Rf_t = Rf_vec(t);                                                  % (1*1) Risk-free rate R^f_t

        date_str = num2str(dates(t));
        R_axis = Smooth_AllR.(date_str);                                   % (1*30000) Gross return grid
        f_star_curve = Smooth_AllR_RND.(date_str);                         % (1*30000) Q-measure PDF


        % === Step 2: build B-spline basis on current grid ===
        min_knot = min(R_axis);
        max_knot = max(R_axis);
        N = numel(R_axis);
        degree = 3;                                                        % cubic B-spline

        B_stack = zeros(b+1, N);
        for i = 1:(b+1)
            B_stack(i, :) = Bspline_basis_function_value( ...
                degree, b, min_knot, max_knot, i, R_axis);
        end

        g_vec = theta * B_stack;                                           % (1×(b+1)) × ((b+1)×N) → 1×N

        integrand = f_star_curve .* g_vec;                                 % (1*30000) f^* × ∑
        integral_val = trapz(R_axis, integrand);                           % (1*1)     ∫ f^* ∑ dR

        delta_t = Rf_t - integral_val;                                     % (1*1)     δ_t = =R^f_t - ∫...
        delta_vec(t) = delta_t;

        denom = delta_t + g_vec;
        minDen = min(denom);
        if minDen <= 0
            warning('t=%d: delta + g_vec has non-positive values (min=%.3e).', t, minDen);
        elseif minDen < 1e-6
            warning('t=%d: delta + g_vec is very close to zero (min=%.3e).', t, minDen);
        end

        % === Step 3: Evaluate M(R_{t+1}; θ) and f_physical curve ===
        M_grid  = 1 ./ (delta_t + g_vec);                                   % (1*30000) Pricing kernel M(R_{t+1}; θ)
        M_vec(t, :) = M_grid;
        f_curve = f_star_curve ./ (Rf_t * M_grid);


        % === Step 4: Interpolate f_t(R_{t+1}; θ) ===
        if R_realized_t < min(R_axis) || R_realized_t > max(R_axis)
            f_physical_at_realized = 1e-10;                                % Penalize out-of-bound realized returns
        else
            f_physical_at_realized = interp1(R_axis, f_curve, R_realized_t, 'linear');
        end


        % === Step 5: Accumulate log-likelihood ===
        LL = LL + log(max(f_physical_at_realized, 1e-12));

    end
end


%% Local Function: non-linear constraint

function [c, ceq] = nlcon_pos(theta, Rf_vec, Smooth_AllR, Smooth_AllR_RND, dates, b, epsPos)

    T = numel(dates);
    c_all = [];

    for t = 1:T
        date_str = num2str(dates(t));
        R_axis   = Smooth_AllR.(date_str);
        fstar    = Smooth_AllR_RND.(date_str);

        degree = 3;  min_knot = min(R_axis);  max_knot = max(R_axis);
        B_stack = zeros(b+1, numel(R_axis));
        for i = 1:(b+1)
            B_stack(i,:) = Bspline_basis_function_value(degree, b, min_knot, max_knot, i, R_axis);
        end
        g_vec = theta * B_stack;

        Rf_t = Rf_vec(t);
        delta_t = Rf_t - trapz(R_axis, fstar .* g_vec);

        denom = delta_t + g_vec;
        minDenom = min(denom);
        c_all = [c_all; epsPos - minDenom];
    end
    c   = c_all;
    ceq = [];
end
