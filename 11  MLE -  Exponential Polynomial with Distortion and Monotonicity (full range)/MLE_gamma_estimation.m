%% Main Function: MLE gamma Estimation for Pricing Kernel

function [gamma_hat, log_lik, delta_vec, M_vec] = MLE_gamma_estimation( ...
    Smooth_AllR, Smooth_AllR_RND, Realized_Return, Risk_Free_Rate, ...
    L, use_delta, alpha, beta, Global_Min_R, Global_Max_R)

    % Settings
    rng(0);
    dates  = Realized_Return.date;
    R_vec  = Realized_Return.realized_ret;
    Rf_vec = Risk_Free_Rate;
    
    % Data Pre-processing
    disp('Pre-processing data for parallel computing...');
    T = length(dates);
    R_grids_All = cell(T, 1);
    f_star_All  = cell(T, 1);
    
    for t = 1:T
        date_str = num2str(dates(t));
        R_grids_All{t} = Smooth_AllR.(date_str);
        f_star_All{t}  = Smooth_AllR_RND.(date_str);
    end

    % Initial guess
    gamma0 = zeros(L,1);
    LB = [];
    UB = [];

    % Optimization options
    options = optimoptions('fmincon', ...
        'Display', 'iter-detailed', ...
        'Algorithm', 'sqp', ...
        'ConstraintTolerance', 1e-9, ...
        'StepTolerance', 1e-9, ...
        'MaxFunctionEvaluations', 50000, ...
        'SpecifyObjectiveGradient', false);

    % Define objective function
    obj_fun = @(param) -log_likelihood_function_par(param, R_vec, Rf_vec, L, ...
        R_grids_All, f_star_All, use_delta, alpha, beta);
    % obj_fun = @(param) -log_likelihood_function(param, R_vec, Rf_vec, L, ...
    %     Smooth_AllR, Smooth_AllR_RND, dates, use_delta, ...
    %     alpha, beta);
    
    % Enforce M_grid decreasing on the FULL Range of R
    nonlcon = @(g) monotone_nonlcon(g, L, Global_Min_R, Global_Max_R);

    % Run optimization
    if isempty(gcp('nocreate')), parpool; end
    disp('Starting optimization with parallel computing...');
    [gamma_hat, neg_LL, exitflag, ~] = fmincon( ...
        obj_fun, gamma0, [], [], [], [], LB, UB, nonlcon, options);

    % Return log-likelihood (positive value)
    log_lik = -neg_LL;

    % Check optimization result
    if     exitflag > 0,  disp('fmincon converged successfully.');
    elseif exitflag == 0, disp('fmincon reached the iteration limit.');
    else,                 disp(['fmincon failed. Exit flag: ', num2str(exitflag)]);
    end
    
    % Display output
    disp('Estimated coefficients:');
    disp(gamma_hat);
    disp(['Final log-likelihood = ', num2str(log_lik)]);

    % Post-estimation call
    [~, delta_vec, M_vec] = log_likelihood_function_par(gamma_hat, R_vec, Rf_vec, L, ...
        R_grids_All, f_star_All, use_delta, alpha, beta);

end