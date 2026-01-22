%% Log-Likelihood Function (Cubic B-Spline Version)

function [LL, delta_vec, M_vec] = log_likelihood_bspline(theta, R_vec, Rf_vec, ...
    Basis_Precomputed, Smooth_AllR, Smooth_AllR_RND, dates, use_delta, alpha, beta)

    T = length(R_vec);
    LL_contributions = zeros(T, 1);
    
    theta = theta(:);
    
    if nargout > 1
        delta_vec = zeros(T, 1);
    else
        delta_vec = [];
    end
    
    if nargout > 2
        date_str_1 = num2str(dates(1));
        temp_grid  = Smooth_AllR.(date_str_1);
        N_grid     = length(temp_grid);
        M_vec      = zeros(T, N_grid);
    else
        M_vec = [];
    end

    % Loop over time
    for t = 1:T
        % === Step 1: Basic inputs ===
        R_realized_t = R_vec(t);
        Rf_t         = Rf_vec(t);
        date_str     = num2str(dates(t));
        
        R_axis       = Smooth_AllR.(date_str);                             % (1 x N)
        f_star_curve = Smooth_AllR_RND.(date_str);                         % (1 x N)
        
        B_mat = Basis_Precomputed{t};                                      % (N x (b+1))
        
        if isempty(B_mat)
             LL_contributions(t) = log(1e-12);
             continue; 
        end
        
        % === Step 2: Compute Spline Sum Q(R) ===
        Spline_Sum = (B_mat * theta)';                                     % (1 x N)
        
        % 數值穩定控制 (避免 exp 爆掉)
        Spline_Sum = max(min(Spline_Sum, 60), -60);
        
        % === Step 3: Compute delta_t (Eq 5) ===
        % Eq (5): delta_t = -ln(Rf) + ln( integral( f* x exp(-Spline_Sum) ) )
        
        if use_delta
            integrand = f_star_curve .* exp(-Spline_Sum);
            integral_val = max(trapz(R_axis, integrand), 1e-300);
            
            delta_t = -log(Rf_t) + log(integral_val);
        else
            delta_t = 0;
        end
        
        if nargout > 1, delta_vec(t) = delta_t; end
        
        % === Step 4: Evaluate M(R_{t+1}) (Eq 4) ===
        % Eq (4): M = exp( delta_t + Spline_Sum )
        
        logM   = delta_t + Spline_Sum;
        M_grid = exp(logM);  
        
        if nargout > 2
            M_vec(t, :) = M_grid;
        end
        
        % === Step 5: Evaluate f_t(R) (Eq 2) ===
        % Eq (2): f_t = f* / (Rf * M)
        % 代入 M: f_t = f* / (Rf * exp(delta + Spline_Sum))
        %             = f* / (Rf * (Integral/Rf) * exp(Spline_Sum))  <-- 因為 exp(delta) = Integral/Rf
        %             = f* / (Integral * exp(Spline_Sum))
        %             = (f* * exp(-Spline_Sum)) / Integral
        
        % 計算 Physical Density (未經 distortion)
        % 這裡直接用上面算好的 integrand (即 f* * exp(-Spline)) 除以 integral_val
        % 這樣數值上最穩
        
        baseline_pdf = integrand ./ integral_val; 
        
        % Regularization check
        if ~all(isfinite(baseline_pdf)) || sum(baseline_pdf)==0
            LL_contributions(t) = log(1e-12);
            continue
        end
        
        % 計算 CDF tildeF
        tildeF = cumtrapz(R_axis, baseline_pdf);
        tildeF = tildeF ./ max(tildeF(end), 1e-12);                        % 確保歸一化
        tildeF = min(max(tildeF, 1e-12), 1-1e-12);                         % 邊界保護
        
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
        % 根據實現的 R_realized_t 插值取出機率密度
        
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
    
    LL = sum(LL_contributions);
end