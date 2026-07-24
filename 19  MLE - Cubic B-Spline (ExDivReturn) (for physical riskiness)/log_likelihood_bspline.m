function [LL, BIC, kappa_vec, SDF_Cell, pit_vec, ...
    LL_contributions, Physical_PDF_Cell] = log_likelihood_bspline( ...
    theta, R_vec, Rf_vec, Basis_Precomputed, R_Grid_Cell, ...
    RND_PDF_Cell, Basis_At_Realized, RND_At_Realized)
%LOG_LIKELIHOOD_BSPLINE No-distortion physical-density likelihood.
%
% q(R)       = B(R) * theta
% kappa_t    = -log(Rf_t) + log integral[f*_t(R) exp(-q(R)) dR]
% M_t(R)     = exp(kappa_t + q(R))
% f_t(R)     = f*_t(R) / (Rf_t M_t(R))
%            = f*_t(R) exp(-q(R)) / integral[f*_t exp(-q) dR]

    theta = theta(:);
    T = numel(R_vec);
    LL_contributions = repmat(log(1e-12), T, 1);

    need_kappa = nargout >= 3;
    need_sdf = nargout >= 4;
    need_pit = nargout >= 5;
    need_physical_pdf = nargout >= 7;

    if need_kappa, kappa_vec = NaN(T, 1); else, kappa_vec = []; end
    if need_sdf, SDF_Cell = cell(T, 1); else, SDF_Cell = {}; end
    if need_pit, pit_vec = NaN(T, 1); else, pit_vec = []; end
    if need_physical_pdf
        Physical_PDF_Cell = cell(T, 1);
    else
        Physical_PDF_Cell = {};
    end

    for t = 1:T
        R_axis = R_Grid_Cell{t};
        rnd_pdf = RND_PDF_Cell{t};
        B_mat = Basis_Precomputed{t};

        spline_sum = B_mat * theta;
        spline_sum = max(min(spline_sum, 60), -60);
        unnormalized_physical = rnd_pdf .* exp(-spline_sum);
        integral_value = trapz(R_axis, unnormalized_physical);

        if ~isfinite(integral_value) || integral_value <= 0
            continue
        end

        kappa_t = -log(Rf_vec(t)) + log(integral_value);
        physical_pdf = unnormalized_physical ./ integral_value;

        if need_kappa
            kappa_vec(t) = kappa_t;
        end
        if need_sdf
            SDF_Cell{t} = exp(kappa_t + spline_sum);
        end
        if need_physical_pdf
            Physical_PDF_Cell{t} = physical_pdf;
        end

        if ~isempty(Basis_At_Realized{t}) && ...
                isfinite(RND_At_Realized(t)) && RND_At_Realized(t) > 0
            spline_at_realized = Basis_At_Realized{t} * theta;
            spline_at_realized = max(min(spline_at_realized, 60), -60);
            density_at_realized = RND_At_Realized(t) * ...
                exp(-spline_at_realized) / integral_value;
            if isfinite(density_at_realized) && density_at_realized > 0
                LL_contributions(t) = log(max(density_at_realized, 1e-12));
            end
        end

        if need_pit
            physical_cdf = cumtrapz(R_axis, physical_pdf);
            if physical_cdf(end) > 0
                physical_cdf = physical_cdf ./ physical_cdf(end);
                pit_value = interp1(R_axis, physical_cdf, R_vec(t), ...
                    'pchip', NaN);
                if isfinite(pit_value)
                    pit_vec(t) = min(max(pit_value, 1e-12), 1 - 1e-12);
                end
            end
        end
    end

    LL = sum(LL_contributions);
    num_free_parameters = numel(theta);
    BIC = num_free_parameters * log(T) - 2 * LL;
end
