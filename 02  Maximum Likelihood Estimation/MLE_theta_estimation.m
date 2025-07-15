%% Main Function: MLE Theta Estimation for Pricing Kernel

function [theta_hat, log_lik] = MLE_theta_estimation(Smooth_AllR, Smooth_AllR_RND, Realized_Return, Risk_Free_Rate, N)

    % Set seed for reproducibility
    rng(0);

    % Extract necessary inputs
    dates   = Realized_Return.date;           % time index
    R_vec   = Realized_Return.realized_ret;   % R_{t+1}
    Rf_vec  = Risk_Free_Rate.rate;            % R^f_t

    % Initial guess for theta (c_1,...,c_N)
    theta0 = zeros(N, 1);

    % Optimization options
    options = optimoptions('fmincon', ...
        'Display', 'iter-detailed', ...
        'Algorithm', 'interior-point', ...
        'SpecifyObjectiveGradient', false);

    % Bounds (optional, currently none)
    LB = [];
    UB = [];

    % Define objective function (fmincon minimizes, so we negate the log-likelihood)
    obj_fun = @(theta) -log_likelihood_function(theta, R_vec, Rf_vec, ...
                        N, Smooth_AllR, Smooth_AllR_RND, dates);

    % Run optimization
    [theta_hat, neg_LL, exitflag, ~] = fmincon(obj_fun, theta0, [], [], [], [], LB, UB, [], options);

    % Return log-likelihood (positive value)
    log_lik = -neg_LL;

    % Check optimization result
    if exitflag > 0
        disp('✅ fmincon converged successfully.');
    elseif exitflag == 0
        disp('⚠️ fmincon reached the iteration limit.');
    else
        disp(['❌ fmincon failed. Exit flag: ', num2str(exitflag)]);
    end
    
    % Display final result
    disp('Estimated parameters (θ):');
    disp(theta_hat);
    disp(['Final log-likelihood = ', num2str(-neg_LL)]);

end


%% Local Function: Log-Likelihood Function

function LL = log_likelihood_function(theta, R_vec, Rf_vec, N, Smooth_AllR, Smooth_AllR_RND, dates)

    T = length(R_vec);
    LL = 0;

    % Parameter vector: theta = [c_1, ..., c_N]
    c = theta(:);

    for t = 1:T

        % === Step 1: Basic inputs ===
        R_realized_t = R_vec(t);                                           % (1*1) Realized gross return R_{t+1}
        Rf_t = Rf_vec(t);                                                  % (1*1) Risk-free rate R^f_t

        date_str = num2str(dates(t));
        R_axis = Smooth_AllR.(date_str);                                   % (1*30000) Gross return grid
        f_star_curve = Smooth_AllR_RND.(date_str);                         % (1*30000) Q-measure PDF
        logR_grid = log(R_axis);                                           % (1*30000) log-return grid for integration


        % === Step 2: Compute δ_t from footnote 4 ===
        poly_sum_neg = zeros(size(logR_grid));                             % (1*30000) -∑_{i=1}^N c_i * (log R)^i
        for i = 1:N
            poly_sum_neg = poly_sum_neg - c(i) * (logR_grid.^i);
        end

        integrand = f_star_curve .* exp(poly_sum_neg);                     % (1*30000) f^* * exp{ -∑ c_i (log R)^i }
        integral_val = trapz(R_axis, integrand);                           % (1*1)     ∫_0^∞ f^* exp{...} dR

        delta_t = -log(Rf_t) + log(integral_val);                          % (1*1)     δ_t = -ln R^f_t + ln ∫...


        % === Step 3: Evaluate M(R_{t+1}; θ) from equation (6) ===
        poly_curve = zeros(size(logR_grid));
        for i = 1:N
            poly_curve = poly_curve + c(i) .* (logR_grid.^i);
        end
        M_grid = exp(delta_t + poly_curve);                                % (1*30000) Pricing kernel M(R_{t+1}; θ)


        % === Step 4: Evaluate f_t(R_{t+1}; θ) from equation (2) ===
        f_physical_curve = f_star_curve ./ (Rf_t * M_grid);                % (1*30000) Conditional physical density


        % === Step 5: Interpolate f_t(R_{t+1}; θ) ===
        if R_realized_t < min(R_axis) || R_realized_t > max(R_axis)
            f_physical_at_realized = 1e-10;                                % Penalize out-of-bound realized returns
        else
            f_physical_at_realized = interp1(R_axis, f_physical_curve, R_realized_t, 'linear');
        end


        % === Step 6: Accumulate log-likelihood ===
        if f_physical_at_realized > 0
            LL = LL + log(f_physical_at_realized);                         % (1*1) log-likelihood
        else
            LL = LL - 1e6;                                                 % Penalize invalid density values
        end

    end
end
