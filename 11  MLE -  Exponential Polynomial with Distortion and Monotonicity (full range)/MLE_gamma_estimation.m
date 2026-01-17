%% Main Function: MLE gamma Estimation for Pricing Kernel

function [gamma_hat, log_lik, delta_vec, M_vec] = MLE_gamma_estimation( ...
    Smooth_AllR, Smooth_AllR_RND, Realized_Return, Risk_Free_Rate, ...
    L, use_delta, alpha, beta, Global_Min_R, Global_Max_R)

    % Settings
    rng(0);
    dates  = Realized_Return.date;
    R_vec  = Realized_Return.realized_ret;
    Rf_vec = Risk_Free_Rate;

    % --- [Step 1] Linear Inequality Constraints ---
    Target_Points = 10002;
    R_check = generate_non_uniform_grid(Global_Min_R, Global_Max_R, Target_Points);
    x_check = log(R_check);
    
    % Basis Derivative Matrix
    % SDF monotone decreasing => P'(ln R) > 0
    % P'(x) = sum( l * gamma_l * x^(l-1) )
    powers = 1:L;
    % Basis_Derivative: l * x^(l-1) each row
    Basis_Derivative = (x_check .^ (powers - 1)) .* powers;
    
    % A*x <= b
    %  Basis_Derivative * gamma >= tol
    % -Basis_Derivative * gamma <= -tol
    tol_strict = 1e-7;
    A_ineq = -Basis_Derivative;
    b_ineq = -tol_strict * ones(size(R_check, 1), 1);

    % --- [Step 2] Optimization Setup ---
    gamma0 = zeros(L,1);
    
    % Optimization options
    options = optimoptions('fmincon', ...
        'Display', 'off', ... 
        'Algorithm', 'sqp', ...
        'ConstraintTolerance', 1e-9, ...
        'StepTolerance', 1e-9, ...
        'MaxFunctionEvaluations', 50000);
        
    % Define objective function
    obj_fun = @(param) -log_likelihood_function(param, R_vec, Rf_vec, L, ...
        Smooth_AllR, Smooth_AllR_RND, dates, use_delta, alpha, beta);

    % --- [Step 3] Run Optimization (Serial) ---
    [gamma_hat, neg_LL, ~] = fmincon( ...
        obj_fun, gamma0, A_ineq, b_ineq, [], [], [], [], [], options);
        
    % Return log-likelihood (positive value)
    log_lik = -neg_LL;
    
    % Post-estimation call
    [~, delta_vec, M_vec] = log_likelihood_function(gamma_hat, R_vec, Rf_vec, L, ...
        Smooth_AllR, Smooth_AllR_RND, dates, use_delta, alpha, beta);

end