%% Main Function: MLE gamma Estimation for Pricing Kernel

function [gamma_hat, log_lik, kappa_vec, M_vec] = MLE_gamma_estimation( ...
    Smooth_AllR, Smooth_AllR_RND, Realized_Return, Risk_Free_Rate, L, use_kappa)

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
        Smooth_AllR, Smooth_AllR_RND, dates, use_kappa);

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
    [~, kappa_vec, M_vec] = log_likelihood_function(gamma_hat, R_vec, Rf_vec, L, ...
    Smooth_AllR, Smooth_AllR_RND, dates, use_kappa);

end


%% Local Function: Log-Likelihood Function

function [LL, kappa_vec, M_vec] = log_likelihood_function(param, R_vec, Rf_vec, L, ...
    Smooth_AllR, Smooth_AllR_RND, dates, use_kappa)

    T = length(R_vec);
    LL = 0;

    % Parameter vector: gamma
    gamma = param(:);

    % Keep kappa_vec & M_vec
    R_axis_1  = Smooth_AllR.(num2str(dates(1)));
    N         = numel(R_axis_1);
    kappa_vec = NaN(T, 1);
    M_vec     = NaN(T, N);

    % log-likelihood
    for t = 1:T

        % === Step 1: Basic inputs ===
        R_realized_t = R_vec(t);                                           % (1*1) Realized gross return R_{t+1}
        Rf_t = Rf_vec(t);                                                  % (1*1) Risk-free rate R^f_t

        date_str = num2str(dates(t));
        R_axis = Smooth_AllR.(date_str);                                   % (1*30000) Gross return grid
        f_star_curve = Smooth_AllR_RND.(date_str);                         % (1*30000) Q-measure PDF        
        logR_grid = log(R_axis);                                           % (1*30000) log-return grid for integration


        % === Step 2: Compute polynomial part ===
        poly_sum = zeros(size(logR_grid));                                 % (1*30000)
        for l = 1:L
            poly_sum = poly_sum + gamma(l) * (logR_grid.^l);
        end
        

        % === Step 3: Compute kappa_t if enabled ===
        if use_kappa
            integrand = f_star_curve .* exp(poly_sum);                     % (1*30000)
            integral_val = trapz(R_axis, integrand);                       % (1*1)
            kappa_t = -log(Rf_t) + log(integral_val);                      % (1*1)
        else
            kappa_t = 0;
        end
        kappa_vec(t) = kappa_t;


        % === Step 4: Evaluate M(R_{t+1}; γ) ===
        logM   = kappa_t - poly_sum;
        M_grid = exp(logM);                                                % (1*30000)
        M_vec(t, :) = M_grid;


        % === Step 5: Evaluate f_t(R; γ) ===
        Rf_t_times_M_grid = Rf_t * M_grid;
        f_physical_curve = f_star_curve ./ Rf_t_times_M_grid;                % (1*30000)

        % Renorm for numerical drift
        Z = trapz(R_axis, f_physical_curve);
        if ~isfinite(Z) || Z<=0, LL = -Inf; return;
        else, f_physical_curve = f_physical_curve ./ Z;
        end


        % === Step 6: Interpolate f_t(R_{t+1}; γ) ===
        if R_realized_t < min(R_axis) || R_realized_t > max(R_axis)
            f_physical_at_realized = 1e-10;                                % Penalize out-of-bound realized returns
        else
            f_physical_at_realized = interp1(R_axis, f_physical_curve, R_realized_t,'pchip', 'extrap');
        end


        % === Step 7: Accumulate log-likelihood ===
        LL = LL + log(max(f_physical_at_realized, 1e-12));

    end
end
