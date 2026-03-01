%% Log-Likelihood Function (Generalized B-Spline Version)

function [LL, BIC, delta_vec, M_vec, pit_vec] = log_likelihood_bspline(theta, R_vec, Rf_vec, ...
    Basis_Precomputed, Smooth_AllR, Smooth_AllR_RND, months, use_delta, alpha, beta)

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
        Spline_Sum = (B_mat * theta);
        Spline_Sum = max(min(Spline_Sum, 60), -60);
        
        % === Step 3: Compute delta_t ===
        % delta_t = -ln(Rf) + ln( integral( f* x exp(-Spline_Sum) ) )
        
        if use_delta
            integrand = f_star_curve .* exp(-Spline_Sum);
            integral_val = max(trapz(R_axis, integrand), 1e-300);
            
            delta_t = -log(Rf_t) + log(integral_val);
        else
            delta_t = 0;
            integrand = f_star_curve .* exp(-Spline_Sum);
            integral_val = 1;
        end
        
        if nargout > 2, delta_vec(t) = delta_t; end
        
        % === Step 4: Evaluate M(R_{t+1}) ===
        % M = exp( delta_t + Spline_Sum )        
        logM   = delta_t + Spline_Sum;
        M_grid = exp(logM);  
        
        if nargout > 3
            M_vec(t, :) = M_grid';
        end
        
        % === Step 5: Evaluate f_t(R) ===
        % f_t = f* / (Rf * M)
        % 代入 M: f_t = f* / (Rf * exp(delta + Spline_Sum))
        %             = f* / (Rf * (Integral/Rf) * exp(Spline_Sum))
        %             = f* / (Integral * exp(Spline_Sum))
        %             = (f* * exp(-Spline_Sum)) / Integral        
        baseline_pdf = integrand ./ integral_val; 
        
        % Regularization check
        if ~all(isfinite(baseline_pdf)) || sum(baseline_pdf)==0
            LL_contributions(t) = log(1e-12);
            continue
        end
        
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
        
        % === Step 6: Distortion (Alpha / Beta) ===
        w    = -log(tildeF);                                               
        Dinv = exp( -(w.^(1/alpha)) / beta );                              
        Jac  = Dinv .* ( w.^(1/alpha - 1) ) ./ ( alpha*beta .* tildeF );   
        
        f_physical_curve = Jac .* baseline_pdf;
        
        % Normalize physical curve again just in case
        Z1 = trapz(R_axis, f_physical_curve);
        if ~isfinite(Z1) || Z1<=0
            LL_contributions(t) = log(1e-12);
            continue
        end
        f_physical_curve = f_physical_curve ./ Z1;
        
        % === Step 7: Interpolate to get Likelihood of Realized Return ===        
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
    
    % 注意：在 P-spline 架構下，受限於 Penalty，真實的有效自由度會小於 length(theta)。
    % 這裡保留長度作為簡化的保守估計，未來模型選擇將主要依賴 Cross-Validation 而非 BIC。
    k   = length(theta);
    n   = T;
    BIC = k * log(n) - 2 * LL;
end