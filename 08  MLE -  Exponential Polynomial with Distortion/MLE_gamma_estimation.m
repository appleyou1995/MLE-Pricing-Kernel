%% Main Function: MLE gamma Estimation for Pricing Kernel

function [gamma_hat, log_lik, delta_vec, M_vec] = MLE_gamma_estimation( ...
    Smooth_AllR, Smooth_AllR_RND, Realized_Return, Risk_Free_Rate, ...
    L, use_delta, alpha, beta)

    % Settings
    rng(0);
    dates  = Realized_Return.date;
    R_vec  = Realized_Return.realized_ret;
    Rf_vec = Risk_Free_Rate;

    % Initial guess
    gamma0 = zeros(L,1);
    LB = [];
    UB = [];

    % Optimization options
    options = optimoptions('fmincon', ...
        'Display', 'iter-detailed', ...
        'Algorithm', 'interior-point', ...
        'SpecifyObjectiveGradient', false);

    % Define objective function (fmincon minimizes, so we negate the log-likelihood)
    obj_fun = @(param) -log_likelihood_function(param, R_vec, Rf_vec, L, ...
        Smooth_AllR, Smooth_AllR_RND, dates, use_delta, ...
        alpha, beta);

    % Run optimization
    [gamma_hat, neg_LL, exitflag, ~] = fmincon(obj_fun, gamma0, [], [], [], [], LB, UB, [], options);

    % Return log-likelihood (positive value)
    log_lik = -neg_LL;

    % Check optimization result
    if     exitflag > 0,  disp('✅ fmincon converged successfully.');
    elseif exitflag == 0, disp('⚠️ fmincon reached the iteration limit.');
    else,                 disp(['❌ fmincon failed. Exit flag: ', num2str(exitflag)]);
    end
    
    % Display output
    disp('Estimated coefficients:');
    disp(gamma_hat);
    disp(['Final log-likelihood = ', num2str(log_lik)]);

    % Post-estimation call
    [~, delta_vec, M_vec] = log_likelihood_function(gamma_hat, R_vec, Rf_vec, L, ...
    Smooth_AllR, Smooth_AllR_RND, dates, use_delta, alpha, beta);

end


%% Local Function: Log-Likelihood Function

function [LL, delta_vec, M_vec] = log_likelihood_function(param, R_vec, Rf_vec, L, ...
    Smooth_AllR, Smooth_AllR_RND, dates, use_delta, alpha, beta)

    T = length(R_vec);
    LL = 0;

    % Parameter vector: gamma
    gamma = param(:);

    % Keep delta_vec & M_vec
    R_axis_1  = Smooth_AllR.(num2str(dates(1)));
    N         = numel(R_axis_1);
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
        logR_grid = log(R_axis);                                           % (1*30000) log-return grid


        % === Step 2: Compute polynomial part ===
        poly_sum = zeros(size(logR_grid));                                 % (1*30000)
        for l = 1:L
            poly_sum = poly_sum + gamma(l) * (logR_grid.^l);
        end
        poly_sum = max(min(poly_sum, 60), -60);
        

        % === Step 3: Compute δ_t if enabled ===
        if use_delta
            integrand = f_star_curve .* exp(poly_sum);                     % (1*30000)
            integral_val = max(trapz(R_axis, integrand), 1e-300);          % (1*1)
            delta_t = -log(Rf_t) + log(integral_val);                      % (1*1)
        else
            delta_t = 0;
        end
        delta_vec(t) = delta_t;


        % === Step 4: Evaluate M(R_{t+1}; γ) ===
        logM   = delta_t - poly_sum;
        M_grid = exp(logM);                                                % (1*30000)
        M_vec(t, :) = M_grid;


        % === Step 5: Evaluate f_t(R; γ) ===
        Rf_t_times_M_grid = Rf_t * M_grid;
        baseline_pdf = f_star_curve ./ Rf_t_times_M_grid;                  % (1*30000)

        % baseline regularization
        if ~all(isfinite(baseline_pdf))
            LL = LL + log(1e-12);
            continue
        end
        Z0 = trapz(R_axis, baseline_pdf);
        if ~isfinite(Z0) || Z0<=0
            LL = LL + log(1e-12);
            continue
        end
        baseline_pdf = baseline_pdf ./ Z0;

        tildeF = cumtrapz(R_axis, baseline_pdf);
        tildeF = tildeF ./ max(tildeF(end), 1e-12);
        tildeF = min(max(tildeF, 1e-12), 1-1e-12);

        % Distortion
        FP_distort_inverse = exp(-((-log(tildeF)).^(1/alpha))/beta);

        % Numerical differentiation + regularization
        f_physical_curve = gradient(FP_distort_inverse, R_axis);
        f_physical_curve = max(f_physical_curve, 0);
        Z1 = trapz(R_axis, f_physical_curve);
        if ~isfinite(Z1) || Z1<=0
            LL = LL + log(1e-12);
            continue
        end
        f_physical_curve = f_physical_curve ./ Z1;

        % check before interpolate
        mask = isfinite(R_axis) & isfinite(f_physical_curve);
        R_axis_good = R_axis(mask);
        fP_good     = f_physical_curve(mask);
        if numel(R_axis_good) < 2
            LL = LL + log(1e-12);
            continue
        end


        % === Step 6: Interpolate f_t(R_{t+1}; γ) ===
        if R_realized_t < R_axis_good(1) || R_realized_t > R_axis_good(end)
            f_physical_at_realized = 1e-12;
        else
            f_physical_at_realized = interp1(R_axis_good, fP_good, R_realized_t, 'pchip');
                if ~isfinite(f_physical_at_realized) || f_physical_at_realized<=0
                    f_physical_at_realized = 1e-12;
                end
        end


        % === Step 7: Accumulate log-likelihood ===
        LL = LL + log(f_physical_at_realized);

    end
end