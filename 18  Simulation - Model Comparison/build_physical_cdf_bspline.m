%% Construct the month-specific physical CDFs implied by a fitted model

function [R_axis_cell, F_physical_cell] = build_physical_cdf_bspline( ...
    theta, Rf_vec, Basis_Precomputed, ...
    Smooth_AllR, Smooth_AllR_RND, months, alpha, beta)

    T = numel(months);
    theta = theta(:);

    if numel(Rf_vec) ~= T || numel(Basis_Precomputed) ~= T
        error('Rf_vec, Basis_Precomputed, and months must have the same length.');
    end

    R_axis_cell     = cell(T, 1);
    F_physical_cell = cell(T, 1);

    for t = 1:T
        col_name = months{t};

        R_axis_t = Smooth_AllR.(col_name);
        fQ_t     = Smooth_AllR_RND.(col_name);
        B_t      = Basis_Precomputed{t};

        R_axis_t = R_axis_t(:);
        fQ_t     = fQ_t(:);

        if isempty(R_axis_t) || isempty(fQ_t) || isempty(B_t)
            error('Missing RND inputs at t=%d (%s).', t, col_name);
        end

        % Same overflow protection as log_likelihood_bspline.m
        spline_sum = B_t * theta;
        spline_sum = max(min(spline_sum, 60), -60);

        integrand = fQ_t .* exp(-spline_sum);
        integral_val = trapz(R_axis_t, integrand);

        if ~isfinite(integral_val) || integral_val <= 0
            error('Invalid normalization integral at t=%d (%s).', ...
                t, col_name);
        end

        % The risk-free rate cancels after substituting delta_t:
        %   Rf^{-1} fQ / M = fQ exp(-spline_sum) / integral_val.
        baseline_pdf = integrand ./ integral_val;

        tildeF = cumtrapz(R_axis_t, baseline_pdf);
        tildeF = tildeF ./ max(tildeF(end), 1e-12);
        tildeF = min(max(tildeF, 1e-12), 1 - 1e-12);

        % Prelec inverse:
        % D^{-1}(x) = exp(-(-log(x))^(1/alpha) / beta).
        F_physical = exp( ...
            -( (-log(tildeF)).^(1 / alpha) ) ./ beta);

        % Enforce a valid numerical CDF for inverse-transform sampling.
        F_physical = cummax(F_physical);
        F_physical(1)   = 0;
        F_physical(end) = 1;

        if any(~isfinite(F_physical))
            error('Non-finite physical CDF at t=%d (%s).', ...
                t, col_name);
        end

        R_axis_cell{t}     = R_axis_t;
        F_physical_cell{t} = F_physical;
    end
end
