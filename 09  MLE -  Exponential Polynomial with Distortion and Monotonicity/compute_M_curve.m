function [R_axis, M_bar] = compute_M_curve(gamma_hat, L, alpha, beta, ...
    R_vec, Rf_vec, dates_vec, Smooth_AllR, Smooth_AllR_RND)

    if istable(R_vec),     R_vec     = table2array(R_vec);     end
    if istable(Rf_vec),    Rf_vec    = table2array(Rf_vec);    end
    if istable(dates_vec), dates_vec = table2array(dates_vec); end

    [~, delta_vec, M_vec] = log_likelihood_function( ...
        gamma_hat, R_vec, Rf_vec, L, ...
        Smooth_AllR, Smooth_AllR_RND, dates_vec, true, alpha, beta); %#ok<ASGLU>

    R_axis = Smooth_AllR.(num2str(dates_vec(1)));
    R_base = R_axis(:)';
    N      = numel(R_base);

    date_fields = arrayfun(@(d) num2str(d), dates_vec, 'UniformOutput', false);
    T = numel(dates_vec);
    M_interp = NaN(T, N);
    for t = 1:T
        R_t = Smooth_AllR.(date_fields{t});
        M_t = M_vec(t, :);
        M_interp(t, :) = interp1(R_t, M_t, R_base, 'pchip', NaN);
    end
    M_bar  = mean(M_interp, 1, 'omitnan');    
end