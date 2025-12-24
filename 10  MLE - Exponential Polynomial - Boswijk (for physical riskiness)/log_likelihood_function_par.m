%% Local Function: Log-Likelihood Function (Parallel)

function [LL, kappa_vec, M_cell] = log_likelihood_function_par(param, R_vec, Rf_vec, L, ...
    R_grids_All, f_star_All)

    T = length(R_vec);
    
    % Parameter vector
    gamma = param(:);
    
    % Pre-allocation
    kappa_vec        = NaN(T, 1);
    M_cell           = cell(T, 1);
    LL_contributions = zeros(T, 1);
    
    % --- Parallel Loop ---
    parfor t = 1:T

        % === Step 1: Basic inputs ===
        R_realized_t = R_vec(t);
        Rf_t         = Rf_vec(t);
        
        R_axis       = R_grids_All{t};      % (1*N_t)
        f_star_curve = f_star_All{t};       % (1*N_t)
        logR_grid    = log(R_axis);

        % === Step 2: Compute polynomial part ===
        poly_sum = zeros(size(logR_grid));
        for l = 1:L
            poly_sum = poly_sum + gamma(l) * (logR_grid.^l);
        end
        
        % === Step 3: Compute kappa_t with Martingale Restriction ===
        % 3.1 Uncorrected Integral
        integrand_raw    = f_star_curve .* exp(poly_sum);
        integral_val_raw = trapz(R_axis, integrand_raw);
        
        % 3.2 Bias Check (E^Q[R] should be Rf)
        EQ_R_biased = trapz(R_axis, f_star_curve .* R_axis);
        
        % 3.3 Correction Ratio (Whiteboard Logic)
        Correction_Ratio = Rf_t / EQ_R_biased;
        
        % 3.4 Linear Scaling
        integral_val = integral_val_raw * Correction_Ratio;
        
        % --- 強制符合理論下界 (Theoretical Floor) ---
        % 當 L=1 且 gamma > 1 時，積分值理論上不能小於 Rf^gamma
        % 只有當 gamma > 1 時才執行此檢查 (避免干擾 gamma < 1 的探索)
        if L == 1 && gamma(1) > 1
             % 理論最小積分值 = Rf ^ gamma
             Theoretical_Min_Integral = Rf_t ^ gamma(1);
             
             % 如果數值積分因為誤差導致小於理論值，強制使用理論值
             % 這等同於把 kappa "Snap" 回 0 (或下界)
             if integral_val < Theoretical_Min_Integral
                 integral_val = Theoretical_Min_Integral;
             end
        end
        % ------------------------------------------------
        
        % 3.5 kappa
        kappa_t = -log(Rf_t) + log(integral_val);
        
        % [Parfor Output]
        kappa_vec(t) = kappa_t;

        % === Step 4: Evaluate M(R_{t+1}; γ) ===
        logM      = kappa_t - poly_sum;
        M_grid    = exp(logM);
        
        % [Parfor Output]
        M_cell{t} = M_grid;

        % === Step 5: Evaluate f_t(R; γ) ===
        % 理論上 M 變大了 Correction_Ratio 倍，f* 是原始的
        % 但下方的 Renormalization (Z) 會自動抵銷這個常數倍數，所以這裡公式不用動
        Rf_t_times_M_grid = Rf_t * M_grid;
        f_physical_curve  = f_star_curve ./ Rf_t_times_M_grid;

        % Renorm for numerical drift (這步很關鍵，確保總機率為 1)
        Z = trapz(R_axis, f_physical_curve);
        
        % 防呆檢查
        if ~isfinite(Z) || Z <= 0
            LL_contributions(t) = -1e10; % Penalty
        else
            f_physical_curve = f_physical_curve ./ Z;
            
            % === Step 6: Interpolate f_t(R_{t+1}; γ) ===
            if R_realized_t < min(R_axis) || R_realized_t > max(R_axis)
                f_physical_at_realized = 1e-10;
            else
                f_physical_at_realized = interp1(R_axis, f_physical_curve, ...
                    R_realized_t, 'pchip', 'extrap');
            end
            
            % === Step 7: Store Log-Likelihood ===
            val = max(f_physical_at_realized, 1e-12);
            LL_contributions(t) = log(val);
        end
    end
    
    % Sum up all contributions
    LL = sum(LL_contributions);
end