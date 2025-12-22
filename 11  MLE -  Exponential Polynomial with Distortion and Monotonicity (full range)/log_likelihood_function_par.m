%% Log-Likelihood Function (Parallel Version)

function [LL, delta_vec, M_vec] = log_likelihood_function_par(param, R_vec, Rf_vec, L, ...
    R_grids_All, f_star_All, use_delta, alpha, beta)
    
    T = length(R_vec);
    
    % Parameter vector: gamma
    gamma = param(:);
    
    % 預先分配空間 (Pre-allocation)
    sample_grid = R_grids_All{1};
    N = numel(sample_grid);
    
    delta_vec = zeros(T, 1);
    M_vec     = zeros(T, N);
    
    % 用來暫存每個時間點的 Log-Likelihood，避免在 parfor 中做 Reduction 競爭
    LL_contributions = zeros(T, 1);
    
    % --- 平行運算 (Parallel Computing) ---
    gamma_const = parallel.pool.Constant(gamma);
    parfor t = 1:T
        % === Step 1: Basic inputs (直接從 Cell 讀取，極快) ===
        R_realized_t = R_vec(t);
        Rf_t         = Rf_vec(t);
        
        R_axis       = R_grids_All{t};
        f_star_curve = f_star_All{t};
        
        logR_grid = log(R_axis);

        % === Step 2: Compute polynomial part (Horner's Method) ===
        local_gamma = gamma_const.Value;
        poly_sum = local_gamma(L) * ones(size(logR_grid));
        for l = L-1:-1:1
            poly_sum = poly_sum .* logR_grid + local_gamma(l);
        end
        poly_sum = poly_sum .* logR_grid;
        poly_sum = max(min(poly_sum, 60), -60);
        
        % === Step 3: Compute δ_t if enabled ===
        if use_delta
            integrand = f_star_curve .* exp(poly_sum);
            integral_val = max(trapz(R_axis, integrand), 1e-300);
            dt = -log(Rf_t) + log(integral_val);
        else
            dt = 0;
        end
        
        % [Parfor Variable] Slicing Assignment for delta_vec
        delta_vec(t) = dt; 

        % === Step 4: Evaluate M(R_{t+1}; γ) ===
        logM   = dt - poly_sum;
        M_grid = exp(logM);
        
        % [Parfor Variable] Slicing Assignment for M_vec
        M_vec(t, :) = M_grid; 

        % === Step 5: Evaluate f_t(R; γ) ===
        Rf_t_times_M_grid = Rf_t * M_grid;
        baseline_pdf = f_star_curve ./ Rf_t_times_M_grid;

        % baseline regularization
        if ~all(isfinite(baseline_pdf))
            LL_contributions(t) = log(1e-12);
            continue;
        end
        
        Z0 = trapz(R_axis, baseline_pdf);
        if ~isfinite(Z0) || Z0 <= 0
            LL_contributions(t) = log(1e-12);
            continue;
        end
        
        baseline_pdf = baseline_pdf ./ Z0;
        tildeF = cumtrapz(R_axis, baseline_pdf);
        tildeF = tildeF ./ max(tildeF(end), 1e-12);
        tildeF = min(max(tildeF, 1e-12), 1-1e-12);

        % Distortion via closed-form Jacobian
        w    = -log(tildeF);
        Dinv = exp( -(w.^(1/alpha)) / beta );
        Jac  = Dinv .* ( w.^(1/alpha - 1) ) ./ ( alpha*beta .* tildeF );
        f_physical_curve = Jac .* baseline_pdf;

        % f_physical_curve regularization
        Z1 = trapz(R_axis, f_physical_curve);
        if ~isfinite(Z1) || Z1 <= 0
            LL_contributions(t) = log(1e-12);
            continue;
        end
        
        f_physical_curve = f_physical_curve ./ Z1;

        % check before interpolate
        mask = isfinite(R_axis) & isfinite(f_physical_curve);
        R_axis_good = R_axis(mask);
        fP_good     = f_physical_curve(mask);

        if numel(R_axis_good) < 2
            LL_contributions(t) = log(1e-12);
            continue;
        end

        % === Step 6: Interpolate f_t(R_{t+1}; γ) ===
        if R_realized_t < R_axis_good(1) || R_realized_t > R_axis_good(end)
            f_physical_at_realized = 1e-12;
        else
            f_physical_at_realized = interp1(R_axis_good, fP_good, R_realized_t, 'pchip');
            if ~isfinite(f_physical_at_realized) || f_physical_at_realized <= 0
                f_physical_at_realized = 1e-12;
            end
        end

        % === Step 7: Store Log-Likelihood ===
        LL_contributions(t) = log(f_physical_at_realized);
    end
    
    LL = sum(LL_contributions);
end