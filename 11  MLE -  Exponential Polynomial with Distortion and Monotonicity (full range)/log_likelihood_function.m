%% Log-Likelihood Function

function [LL, BIC, delta_vec, M_vec, pit_vec] = log_likelihood_function(param, R_vec, Rf_vec, L, ...
    Smooth_AllR, Smooth_AllR_RND, dates, use_delta, alpha, beta)

    T = length(R_vec);
    LL_contributions = zeros(T, 1);

    % Parameter vector: gamma
    gamma = param(:);

    % Keep delta_vec
    if nargout > 2
        delta_vec = zeros(T, 1);
    else
        delta_vec = [];
    end
    
    % Keep M_vec
    if nargout > 3
        date_str_1 = num2str(dates(1));
        temp_grid  = Smooth_AllR.(date_str_1);
        N_grid     = length(temp_grid);
        
        M_vec = zeros(T, N_grid);
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
        R_realized_t = R_vec(t);                                           % (1*1) Realized gross return R_{t+1}
        Rf_t         = Rf_vec(t);                                          % (1*1) Risk-free rate R^f_t

        date_str     = num2str(dates(t));
        R_axis       = Smooth_AllR.(date_str);                             % (1*30000) Gross return grid
        f_star_curve = Smooth_AllR_RND.(date_str);                         % (1*30000) Q-measure PDF        
        logR_grid    = log(R_axis);                                        % (1*30000) log-return grid


        % === Step 2: Compute polynomial part (Horner's Method) ===
        % P(x) = g1*x + g2*x^2 + ... = x * (g1 + x * (g2 + ...))
        poly_sum = gamma(L) * ones(size(logR_grid));
        for l = L-1:-1:1
            poly_sum = poly_sum .* logR_grid + gamma(l);
        end
        poly_sum = poly_sum .* logR_grid;
        poly_sum = max(min(poly_sum, 60), -60);
        

        % === Step 3: Compute δ_t ===
        if use_delta
            integrand = f_star_curve .* exp(poly_sum);
            integral_val_raw = max(trapz(R_axis, integrand), 1e-300);
            
            % Bias Correction Logic
            EQ_R_biased = trapz(R_axis, f_star_curve .* R_axis);
            Correction_Ratio = Rf_t / EQ_R_biased;
            
            integral_val = integral_val_raw * Correction_Ratio;
            
            % Theoretical Floor
            if L == 1 && gamma(1) > 1
                 Theoretical_Min = Rf_t ^ gamma(1);
                 if integral_val < Theoretical_Min
                     integral_val = Theoretical_Min;
                 end
            end
            
            delta_t = -log(Rf_t) + log(integral_val);
        else
            delta_t = 0;
        end
        
        if nargout > 2, delta_vec(t) = delta_t; end


        % === Step 4: Evaluate M(R_{t+1}; γ) ===
        logM   = delta_t - poly_sum;
        M_grid = exp(logM);                                                % (1*30000)
        if nargout > 3
            M_vec(t, :) = M_grid;
        end

        % === Step 5: Evaluate f_t(R; γ) ===
        Rf_t_times_M_grid = Rf_t * M_grid;
        baseline_pdf = f_star_curve ./ Rf_t_times_M_grid;                  % (1*30000)

        % Regularization
        if ~all(isfinite(baseline_pdf))
            LL_contributions(t) = log(1e-12);
            continue
        end
        Z0 = trapz(R_axis, baseline_pdf);
        if ~isfinite(Z0) || Z0<=0
            LL_contributions(t) = log(1e-12);
            continue
        end
        baseline_pdf = baseline_pdf ./ Z0;

        tildeF = cumtrapz(R_axis, baseline_pdf);
        tildeF = tildeF ./ max(tildeF(end), 1e-12);
        tildeF = min(max(tildeF, 1e-12), 1-1e-12);

        % Physical Probability Integral Transform
        if nargout > 4
            u_tilde = interp1(R_axis, tildeF, R_realized_t, 'pchip');            
            u_tilde = min(max(u_tilde, 1e-12), 1-1e-12);
            w_val   = -log(u_tilde);
            pit_val = exp( -(w_val^(1/alpha)) / beta );            
            pit_vec(t) = pit_val;
        end

        % Distortion via closed-form Jacobian
        w    = -log(tildeF);                                               % w = -ln h(R)
        Dinv = exp( -(w.^(1/alpha)) / beta );                              % D^{-1}(h)
        Jac  = Dinv .* ( w.^(1/alpha - 1) ) ./ ( alpha*beta .* tildeF );   % d D^{-1}/d h

        f_physical_curve = Jac .* baseline_pdf;

        % f_physical_curve regularization
        Z1 = trapz(R_axis, f_physical_curve);
        if ~isfinite(Z1) || Z1<=0
            LL_contributions(t) = log(1e-12);
            continue
        end
        f_physical_curve = f_physical_curve ./ Z1;

        % === Step 6: Interpolate ===
        mask = isfinite(R_axis) & isfinite(f_physical_curve);
        R_axis_good = R_axis(mask);
        fP_good     = f_physical_curve(mask);
        
        if numel(R_axis_good) < 2
            val = 1e-12;
        else
            if R_realized_t < R_axis_good(1) || R_realized_t > R_axis_good(end)
                val = 1e-12;
            else
                val = interp1(R_axis_good, fP_good, R_realized_t, 'pchip');
                if ~isfinite(val) || val<=0, val = 1e-12; end
            end
        end
        LL_contributions(t) = log(val);
    end
    
    LL  = sum(LL_contributions);
    k   = length(gamma);
    n   = T;
    BIC = k * log(n) - 2 * LL;
end