function validation_loss = Compute_validation_error( ...
    gamma_hat, L, Smooth_AllR, Smooth_AllR_RND, ...
    Realized_Return_valid, Risk_Free_Rate_valid, use_delta, ...
    alpha, beta)

    % input
    dates  = Realized_Return_valid.date;
    R_vec  = Realized_Return_valid.realized_ret;
    Rf_vec = Risk_Free_Rate_valid;

    T2   = numel(R_vec);
    Uvec = nan(T2,1);

    % Target uniform moments
    k_list = [3 4 5 6];
    m_tar  = 1 ./ (k_list + 1);

    for t = 1:T2

        % --- Basic inputs for current t ---
        R_realized_t = R_vec(t);
        Rf_t         = Rf_vec(t);

        date_str     = num2str(dates(t));
        R_axis       = Smooth_AllR.(date_str);
        f_star_curve = Smooth_AllR_RND.(date_str);
        logR_grid    = log(R_axis);

        % --- polynomial part ---
        poly_sum = zeros(size(logR_grid));
        for l = 1:L
            poly_sum = poly_sum + gamma_hat(l) * (logR_grid.^l);
        end
        poly_sum = max(min(poly_sum, 60), -60);

        % --- δ_t ---
        if use_delta
            integrand    = f_star_curve .* exp(poly_sum);
            integral_val = max(trapz(R_axis, integrand), 1e-300);
            delta_t      = -log(Rf_t) + log(integral_val);
        else
            delta_t      = 0;
        end

        % --- M(R;γ) ---
        logM   = delta_t - poly_sum;
        M_grid = exp(logM);

        % --- baseline pdf & CDF h(R) ---
        baseline_pdf = f_star_curve ./ (Rf_t * M_grid);
        if ~all(isfinite(baseline_pdf)),  continue; end
        Z0 = trapz(R_axis, baseline_pdf);
        if ~isfinite(Z0) || Z0<=0,        continue; end
        baseline_pdf = baseline_pdf ./ Z0;

        h_grid = cumtrapz(R_axis, baseline_pdf);
        h_grid = h_grid ./ max(h_grid(end), 1e-12);
        h_grid = min(max(h_grid, 1e-12), 1-1e-12);

        if R_realized_t < R_axis(1) || R_realized_t > R_axis(end), continue; end

        % --- U_t = D^{-1}( h_t(R_{t+1}) ) ---
        h_at_R = interp1(R_axis, h_grid, R_realized_t, 'pchip');
        if ~isfinite(h_at_R), continue; end
        h_at_R = min(max(h_at_R, 1e-12), 1-1e-12);

        % Prelec inverse distortion
        Uvec(t) = exp( - ((-log(h_at_R))^(1/alpha)) / beta );
    end

    % delete NaN
    Uvec = Uvec(isfinite(Uvec));
    if isempty(Uvec)
        validation_loss = inf; 
        return
    end

    % score
    m_emp = arrayfun(@(k) mean(Uvec.^k), k_list);
    validation_loss = sum((m_emp - m_tar).^2);
end
