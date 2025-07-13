function LL = log_likelihood_function(theta, R_vec, Rf_vec, f_star_interp, N, Smooth_AllR, Smooth_AllR_RND, dates)

    T = length(R_vec);
    LL = 0;

    % Parameter vector: theta = [c_1, ..., c_N]
    c = theta(:);

    for t = 1:T

        % === Step 1: Basic inputs ===
        R_t = R_vec(t);                                                    % (1*1) Realized gross return R_{t+1}
        logR_t = log(R_t);                                                 % (1*1) log(R_{t+1})
        Rf_t = Rf_vec(t);                                                  % (1*1) Risk-free rate R^f_t
        f_star_t = f_star_interp(t);                                       % (1*1) f^*_t(R_{t+1})

        date_str = num2str(dates(t));
        R_axis = Smooth_AllR.(date_str);                                   % (1*30000) Gross return grid
        f_star_curve = Smooth_AllR_RND.(date_str);                         % (1*30000) Q-measure PDF
        logR_grid = log(R_axis);                                           % (1*30000) log-return grid for integration


        % === Step 2: Compute δ_t ===
        poly_sum_neg = zeros(size(logR_grid));                             % (1*30000) -∑_{i=1}^N c_i * (log R)^i
        for i = 1:N
            poly_sum_neg = poly_sum_neg - c(i) * (logR_grid.^i);
        end

        integrand = f_star_curve .* exp(poly_sum_neg);                     % (1*30000) f^* * exp{ -∑ c_i (log R)^i }
        integral_val = trapz(R_axis, integrand);                           % (1*1)     ∫_0^∞ f^* exp{...} dR

        delta_t = -log(Rf_t) + log(integral_val);                          % (1*1)     δ_t = -ln R^f_t + ln ∫...


        % === Step 3: Evaluate M(R_{t+1}; θ) ===
        poly_val = 0;
        for i = 1:N
            poly_val = poly_val + c(i) * logR_t^i;
        end
        M_val = exp(delta_t + poly_val);                                   % (1*1) Pricing kernel M(R_{t+1}; θ)


        % === Step 4: Evaluate f_t(R_{t+1}; θ) from equation (2) ===
        f_val = f_star_t / (Rf_t * M_val);                                 % (1*1)


        % === Step 5: Accumulate log-likelihood ===
        if f_val > 0
            LL = LL + log(f_val);                                          % (1*1) log-likelihood
        else
            LL = LL - 1e6;                                                 % Penalize invalid density values
        end

    end
end
