%% Main Function: MLE Theta Estimation for Pricing Kernel

function [gamma_hat, log_lik, delta_vec, M_vec] = MLE_theta_estimation( ...
    Smooth_AllR, Smooth_AllR_RND, Realized_Return, Risk_Free_Rate, L, use_delta)

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
        Smooth_AllR, Smooth_AllR_RND, dates, use_delta);

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
    Smooth_AllR, Smooth_AllR_RND, dates, use_delta);

end


%% Local Function: Log-Likelihood Function

function [LL, delta_vec, M_vec] = log_likelihood_function(param, R_vec, Rf_vec, L, ...
    Smooth_AllR, Smooth_AllR_RND, dates, use_delta)

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
        
        % R_axis(~isfinite(R_axis)) = eps;
        % R_axis = max(R_axis, eps);
        % f_star_curve(~isfinite(f_star_curve)) = 0;
        % f_star_curve = max(f_star_curve, 0);
        
        logR_grid = log(R_axis);                                           % (1*30000) log-return grid for integration


        % Check -----------------------------------------------------------
        % if any(~isfinite(logR_grid))
        %     fprintf('[WARN][t=%d %s] logR_grid 含 NaN/Inf；min=%.3g, max=%.3g\n', ...
        %         t, date_str, min(logR_grid), max(logR_grid));
        % end
        % logR_rng = max(logR_grid) - min(logR_grid);
        % if logR_rng > 50
        %     fprintf('[WARN][t=%d %s] logR 範圍極大：range=%.3g（請檢查 R_axis 邊界/單位）\n', ...
        %         t, date_str, logR_rng);
        % end
        % -----------------------------------------------------------------


        % === Step 2: Compute polynomial part ===
        poly_sum = zeros(size(logR_grid));                                 % (1*30000)
        for l = 1:L
            poly_sum = poly_sum + gamma(l) * (logR_grid.^l);
        end
        
        % Check -----------------------------------------------------------
        if any(~isfinite(poly_sum))
            fprintf('[WARN][t=%d %s] poly_sum 出現 NaN/Inf\n', t, date_str);
        end
        % -----------------------------------------------------------------
        

        % === Step 3: Compute δ_t if enabled ===
        if use_delta
            integrand = f_star_curve .* exp(poly_sum);                     % (1*30000)
            integral_val = trapz(R_axis, integrand);                       % (1*1)
            delta_t = -log(Rf_t) + log(integral_val);                      % (1*1)
        else
            delta_t = 0;
        end
        delta_vec(t) = delta_t;

        % === Step 4: Evaluate M(R_{t+1}; γ) ===
        logM   = delta_t - poly_sum;
        % logM   = min(max(logM, -700), 700);
        M_grid = exp(logM);                                                % (1*30000)
        M_vec(t, :) = M_grid;


        % === Step 5: Evaluate f_t(R; γ) ===
        Rf_t_times_M_grid = Rf_t * M_grid;
        % denom = max(Rf_t_times_M_grid, 1e-300);
        f_physical_curve = f_star_curve ./ Rf_t_times_M_grid;                % (1*30000)
        % f_physical_curve(~isfinite(f_physical_curve)) = 0;

        % Z = trapz(R_axis, f_physical_curve);
        % if ~isfinite(Z) || Z <= 0
        %     LL = LL + log(1e-300);
        %     continue
        % end

        % Renorm for numerical drift
        Z = trapz(R_axis, f_physical_curve);
        if ~isfinite(Z) || Z<=0, LL = -Inf; return;
        else, f_physical_curve = f_physical_curve ./ Z;
        end

        % if numel(R_axis) < 2 || numel(f_physical_curve) < 2
        %     fprintf('[ERROR][t=%d] R_axis 或 f_physical_curve 少於2個點\n', t);
        %     keyboard;  % 暫停在這裡
        % end


        % === Step 6: Interpolate f_t(R_{t+1}; γ) ===
        if R_realized_t < min(R_axis) || R_realized_t > max(R_axis)
            f_physical_at_realized = 1e-10;                                % Penalize out-of-bound realized returns
        else
            f_physical_at_realized = interp1(R_axis, f_physical_curve, R_realized_t,'pchip', 'extrap');
        end


        % === Step 7: Accumulate log-likelihood ===
        % val = f_physical_at_realized;
        % val(~isfinite(val) | val<=0) = 1e-300;
        % LL  = LL + log(min(val, 1e300));
        LL = LL + log(max(f_physical_at_realized, 1e-12));

    end
end
