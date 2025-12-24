%% Main Function: MLE gamma Estimation for Pricing Kernel

function [gamma_hat, log_lik, kappa_vec, M_cell] = MLE_gamma_estimation( ...
    Smooth_AllR, Smooth_AllR_RND, Realized_Return, Risk_Free_Rate, L)

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
        Smooth_AllR, Smooth_AllR_RND, dates);

    % Run optimization
    [gamma_hat, neg_LL, exitflag, ~] = fmincon(obj_fun, gamma0, [], [], [], [], LB, UB, [], options);

    % Return log-likelihood (positive value)
    log_lik = -neg_LL;

    % Check optimization result
    if     exitflag > 0,  disp('fmincon converged successfully.');
    elseif exitflag == 0, disp('fmincon reached the iteration limit.');
    else,                 disp(['fmincon failed. Exit flag: ', num2str(exitflag)]);
    end
    
    % Display output
    disp('Estimated coefficients:');
    disp(gamma_hat);
    disp(['Final log-likelihood = ', num2str(log_lik)]);

    % Post-estimation call
    [~, kappa_vec, M_cell] = log_likelihood_function(gamma_hat, R_vec, Rf_vec, L, ...
    Smooth_AllR, Smooth_AllR_RND, dates);

end


%% Local Function: Log-Likelihood Function

function [LL, kappa_vec, M_cell] = log_likelihood_function(param, R_vec, Rf_vec, L, ...
    Smooth_AllR, Smooth_AllR_RND, dates)

    T = length(R_vec);
    LL = 0;

    % Parameter vector: gamma
    gamma = param(:);

    % Keep kappa_vec & M_vec
    % R_axis_1  = Smooth_AllR.(num2str(dates(1)));
    % N         = numel(R_axis_1);
    kappa_vec = NaN(T, 1);
    M_cell    = cell(T, 1);

    % log-likelihood
    for t = 1:T

        % === Step 1: Basic inputs ===
        R_realized_t = R_vec(t);                                           % (1*1) Realized gross return R_{t+1}
        Rf_t         = Rf_vec(t);                                          % (1*1) Risk-free rate R^f_t

        date_str     = num2str(dates(t));
        R_axis       = Smooth_AllR.(date_str);                             % (1*N_t) Gross return grid
        f_star_curve = Smooth_AllR_RND.(date_str);                         % (1*N_t) Q-measure PDF        
        logR_grid    = log(R_axis);                                        % (1*N_t) log-return grid for integration


        % === Step 2: Compute polynomial part ===
        poly_sum = zeros(size(logR_grid));                                 % (1*N_t)
        for l = 1:L
            poly_sum = poly_sum + gamma(l) * (logR_grid.^l);
        end
        

        % === Step 3: Compute kappa_t if enabled ===

        % 3.1 Uncorrected Integral
        integrand_raw    = f_star_curve .* exp(poly_sum);                  % (1*N_t)
        integral_val_raw = trapz(R_axis, integrand_raw);                   % (1*1)

        % 3.2 Bias Check
        EQ_R_biased = trapz(R_axis, f_star_curve .* R_axis);

        % 3.3 Correction Ratio
        % Refined_Density = Ratio * Original_Density
        Correction_Ratio = Rf_t / EQ_R_biased;

        % 3.4 Linear Scaling
        integral_val = integral_val_raw * Correction_Ratio;
        
        % 3.5 kappa
        kappa_t      = -log(Rf_t) + log(integral_val);                     % (1*1)
        kappa_vec(t) = kappa_t;


        % === Step 4: Evaluate M(R_{t+1}; γ) ===
        logM      = kappa_t - poly_sum;
        M_grid    = exp(logM);                                             % (1*N_t)
        M_cell{t} = M_grid;


        % === Step 5: Evaluate f_t(R; γ) ===
        Rf_t_times_M_grid = Rf_t * M_grid;
        f_physical_curve  = f_star_curve ./ Rf_t_times_M_grid;             % (1*N_t)

        % Renorm for numerical drift
        Z = trapz(R_axis, f_physical_curve);
        if ~isfinite(Z) || Z<=0, LL = -Inf; return;
        else, f_physical_curve = f_physical_curve ./ Z;
        end


        % === Step 6: Interpolate f_t(R_{t+1}; γ) ===
        if R_realized_t < min(R_axis) || R_realized_t > max(R_axis)
            % Penalize out-of-bound realized returns
            f_physical_at_realized = 1e-10;
        else
            f_physical_at_realized = interp1(R_axis, f_physical_curve, ...
                R_realized_t,'pchip', 'extrap');
        end


        % === Step 7: Accumulate log-likelihood ===
        LL = LL + log(max(f_physical_at_realized, 1e-12));

    end
end
