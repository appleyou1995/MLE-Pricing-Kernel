%% Log-Likelihood Function

function [LL, BIC, delta_vec, M_vec, pit_vec] = log_likelihood_bspline(theta, R_vec, Rf_vec, ...
    Basis_Precomputed, Smooth_AllR, Smooth_AllR_RND, months, alpha, beta)

    T = length(R_vec);
    LL_contributions = zeros(T, 1);
    
    theta = theta(:);
    
    % Keep delta_vec
    if nargout > 2
        delta_vec = zeros(T, 1);
    else
        delta_vec = [];
    end
    
    % Keep M_vec
    if nargout > 3
        col_name_1 = months{1};
        temp_grid  = Smooth_AllR.(col_name_1);
        N_grid     = length(temp_grid);
        M_vec      = zeros(T, N_grid);
    else
        M_vec = [];
    end

    % Keep pit_vec
    if nargout > 4
        pit_vec = zeros(T, 1);
    else
        pit_vec = [];
    end

    % Loop over time
    for t = 1:T
        % === Step 1: Basic inputs ===
        R_realized_t = R_vec(t);
        Rf_t         = Rf_vec(t);
        col_name     = months{t};
        
        R_axis = Smooth_AllR.(col_name);
        f_star_curve = Smooth_AllR_RND.(col_name);
        
        B_mat = Basis_Precomputed{t};
        
        if isempty(R_axis) || isempty(f_star_curve) || isempty(B_mat)
             LL_contributions(t) = log(1e-12);
             continue; 
        end

        R_axis = R_axis(:);
        f_star_curve = f_star_curve(:);
        
        % === Step 2: Compute Spline Sum Q(R) ===
        Spline_Sum = B_mat * theta;

        % === Step 3: Compute delta_t in the log domain ===
        % log_integrand = log(f*) - Spline_Sum
        % log_integral  = log integral exp(log_integrand) dR
        if any(f_star_curve < 0) || any(~isfinite(f_star_curve))
            LL_contributions(t) = log(1e-12);
            continue
        end

        log_f_star = log(f_star_curve);
        log_integrand = log_f_star - Spline_Sum;
        log_integral = log_trapz_exp(R_axis, log_integrand);

        if ~isfinite(log_integral)
            LL_contributions(t) = log(1e-12);
            continue
        end

        delta_t = -log(Rf_t) + log_integral;

        if nargout > 2
            delta_vec(t) = delta_t;
        end

        % === Step 4: Evaluate M(R) ===
        logM = delta_t + Spline_Sum;

        if nargout > 3
            M_grid = exp(logM);
            M_vec(t, :) = M_grid';
        end

        % === Step 5: Evaluate baseline physical density ===
        log_baseline_pdf = log_integrand - log_integral;
        baseline_pdf = exp(log_baseline_pdf);

        % Keep the existing numerical guard unchanged
        if ~all(isfinite(baseline_pdf)) || sum(baseline_pdf) == 0
            LL_contributions(t) = log(1e-12);
            continue
        end

        tildeF = cumtrapz(R_axis, baseline_pdf);
        tildeF = tildeF ./ max(tildeF(end), 1e-12);
        tildeF = min(max(tildeF, 1e-12), 1 - 1e-12);

        % Physical Probability Integral Transform
        if nargout > 4
            u_tilde = interp1( ...
                R_axis, tildeF, R_realized_t, 'pchip');

            u_tilde = min(max(u_tilde, 1e-12), 1 - 1e-12);
            w_val = -log(u_tilde);

            pit_val = exp( ...
                -(w_val^(1 / alpha)) / beta);

            pit_vec(t) = pit_val;
        end

        % === Step 6: Distortion in the log domain ===
        w = -log(tildeF);

        log_Dinv = -(w.^(1 / alpha)) ./ beta;

        log_Jac = log_Dinv + ...
            (1 / alpha - 1) .* log(w) - ...
            log(alpha) - log(beta) - log(tildeF);

        log_f_physical_unnormalized = ...
            log_Jac + log_baseline_pdf;

        log_Z1 = log_trapz_exp( ...
            R_axis, log_f_physical_unnormalized);

        if ~isfinite(log_Z1)
            LL_contributions(t) = log(1e-12);
            continue
        end

        log_f_physical_curve = ...
            log_f_physical_unnormalized - log_Z1;

        % === Step 7: Interpolate log density at realized return ===
        mask = isfinite(R_axis) & ...
            isfinite(log_f_physical_curve);

        R_axis_good = R_axis(mask);
        log_fP_good = log_f_physical_curve(mask);

        if numel(R_axis_good) < 2
            log_val = log(1e-12);
        else
            if R_realized_t < R_axis_good(1) || ...
                    R_realized_t > R_axis_good(end)

                log_val = log(1e-12);
            else
                log_val = interp1( ...
                    R_axis_good, ...
                    log_fP_good, ...
                    R_realized_t, ...
                    'pchip');

                if ~isfinite(log_val)
                    log_val = log(1e-12);
                end
            end
        end

        LL_contributions(t) = log_val;
    end
    
    LL  = sum(LL_contributions);
    k   = length(theta);
    n   = T;
    BIC = k * log(n) - 2 * LL;
end