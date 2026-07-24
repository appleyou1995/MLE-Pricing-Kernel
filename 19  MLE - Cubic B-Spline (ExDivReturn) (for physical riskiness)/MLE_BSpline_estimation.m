function [theta_hat, log_lik, BIC, exitflag, output, kappa_vec, ...
    SDF_Cell, pit_vec, LL_contributions, Model_Info] = ...
    MLE_BSpline_estimation(Smooth_AllR, Smooth_AllR_RND, ...
    Realized_Return, Risk_Free_Rate, Num_Basis_Function, ...
    Spline_Degree, Enforce_Decreasing, Global_Min_R, Global_Max_R)
%MLE_BSPLINE_ESTIMATION Estimate a cubic B-spline SDF without distortion.
%
% log M_t(R) = kappa_t + B(R) * theta

    rng(0);
    R_vec = double(Realized_Return.realized_ret(:));
    Rf_vec = double(Risk_Free_Rate(:));
    months = Smooth_AllR.Properties.VariableNames;
    rnd_months = Smooth_AllR_RND.Properties.VariableNames;
    T = numel(R_vec);

    if numel(Rf_vec) ~= T || numel(months) ~= T || numel(rnd_months) ~= T
        error('Input lengths are inconsistent in MLE_BSpline_estimation.');
    end
    if Spline_Degree ~= 3
        error('This program is intentionally restricted to cubic splines (degree 3).');
    end

    Spline_Order = Spline_Degree + 1;
    Num_Breaks = Num_Basis_Function - Spline_Order + 2;
    if Num_Breaks < 2
        error(['Num_Basis_Function = %d is too small for spline degree %d. ' ...
            'At least %d basis functions are required.'], ...
            Num_Basis_Function, Spline_Degree, Spline_Order);
    end

    Breaks = linspace(Global_Min_R, Global_Max_R, Num_Breaks);
    Knots = augknt(Breaks, Spline_Order);

    Basis_Precomputed = cell(T, 1);
    R_Grid_Cell = cell(T, 1);
    RND_PDF_Cell = cell(T, 1);
    Basis_At_Realized = cell(T, 1);
    RND_At_Realized = NaN(T, 1);

    for t = 1:T
        R_axis = double(Smooth_AllR.(months{t})(:));
        rnd_pdf = double(Smooth_AllR_RND.(rnd_months{t})(:));

        R_Grid_Cell{t} = R_axis;
        RND_PDF_Cell{t} = rnd_pdf;
        Basis_Precomputed{t} = spcol(Knots, Spline_Order, R_axis);

        if R_vec(t) >= R_axis(1) && R_vec(t) <= R_axis(end)
            Basis_At_Realized{t} = spcol(Knots, Spline_Order, R_vec(t));
            RND_At_Realized(t) = interp1( ...
                R_axis, rnd_pdf, R_vec(t), 'pchip');
        else
            Basis_At_Realized{t} = [];
        end
    end

    theta0 = zeros(Num_Basis_Function, 1);

    if Enforce_Decreasing
        % theta(i+1) - theta(i) <= 0 is sufficient for a non-increasing
        % B-spline curve and therefore a non-increasing SDF in R.
        A_ineq = zeros(Num_Basis_Function - 1, Num_Basis_Function);
        for i = 1:(Num_Basis_Function - 1)
            A_ineq(i, i) = -1;
            A_ineq(i, i + 1) = 1;
        end
        b_ineq = zeros(Num_Basis_Function - 1, 1);
    else
        A_ineq = [];
        b_ineq = [];
    end

    options = optimoptions('fmincon', ...
        'Display', 'iter', ...
        'Algorithm', 'sqp', ...
        'ConstraintTolerance', 1e-8, ...
        'OptimalityTolerance', 1e-7, ...
        'StepTolerance', 1e-9, ...
        'MaxIterations', 2000, ...
        'MaxFunctionEvaluations', 50000);

    objective = @(theta) -log_likelihood_bspline( ...
        theta, R_vec, Rf_vec, Basis_Precomputed, R_Grid_Cell, ...
        RND_PDF_Cell, Basis_At_Realized, RND_At_Realized);

    try
        [theta_hat, negative_log_lik, exitflag, output] = fmincon( ...
            objective, theta0, A_ineq, b_ineq, [], [], ...
            [], [], [], options);
    catch ME
        error('fmincon failed before producing a usable estimate: %s', ME.message);
    end

    log_lik = -negative_log_lik;
    [~, BIC, kappa_vec, SDF_Cell, pit_vec, LL_contributions] = ...
        log_likelihood_bspline(theta_hat, R_vec, Rf_vec, ...
            Basis_Precomputed, R_Grid_Cell, RND_PDF_Cell, ...
            Basis_At_Realized, RND_At_Realized);

    Model_Info = struct();
    Model_Info.spline_degree = Spline_Degree;
    Model_Info.spline_order = Spline_Order;
    Model_Info.num_basis_functions = Num_Basis_Function;
    Model_Info.num_free_theta_parameters = Num_Basis_Function;
    Model_Info.breaks = Breaks(:);
    Model_Info.knots = Knots(:);
    Model_Info.global_min_R = Global_Min_R;
    Model_Info.global_max_R = Global_Max_R;
    Model_Info.enforce_decreasing = logical(Enforce_Decreasing);
    Model_Info.distortion_applied = false;
    Model_Info.theta_normalization = "none";
end
