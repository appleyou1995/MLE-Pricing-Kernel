%% Main Function: MLE Theta Estimation for B-Spline Pricing Kernel

function [theta_hat, log_lik, delta_vec, M_vec] = MLE_BSpline_estimation( ...
    Smooth_AllR, Smooth_AllR_RND, Realized_Return, Risk_Free_Rate, ...
    b, use_delta, alpha, beta, Global_Min_R, Global_Max_R)

    % Settings
    rng(0);
    dates  = Realized_Return.date;
    R_vec  = Realized_Return.realized_ret;
    Rf_vec = Risk_Free_Rate;
    months = Smooth_AllR.Properties.VariableNames;
    T      = length(R_vec);

    % --- [Step 1] Construct Knots & Precompute Basis ---
    n         = 3;
    k_order   = n + 1;
    num_knots = n + b + 2;
    
    min_knot = Global_Min_R;
    max_knot = Global_Max_R;
    
    knots = linspace(min_knot, max_knot, num_knots);
    knots(1:(n+1))      = min_knot; 
    knots((end-n):end)  = max_knot;
    
    % Precompute B-Spline Basis for each day to speed up fmincon
    Basis_Precomputed = cell(T, 1);
    
    for t = 1:T
        try
            col_name = months{t};
            R_axis = Smooth_AllR.(col_name);
            
            % 確保為行向量
            R_axis = R_axis(:);
            
            % 計算 B-spline matrix
            B = spcol(knots, k_order, R_axis);
            
            Basis_Precomputed{t} = B;
            
        catch ME
            warning('Error in precompute at t=%d (%s): %s', t, months{t}, ME.message);
            Basis_Precomputed{t} = [];
        end
    end

    % --- [Step 2] Inequality Constraints (Monotonicity) ---
    % 目標：Pricing Kernel 隨 R 遞減 => exp(delta + sum theta*B) 遞減
    % 因為 B-spline basis 是局部支撐且有序，若 theta_i >= theta_{i+1}，則曲線大致遞減。
    % 建立差分矩陣 D，使得 D * theta <= 0 代表 theta_i - theta_{i-1} <= 0 (遞減)
    % 或是 theta_{i+1} - theta_i <= 0
    
    num_params = b + 1;
    % D * theta <= 0
    % Row i: theta_{i+1} - theta_i <= 0
    D = zeros(b, num_params);
    for i = 1:b
        D(i, i)   = -1;
        D(i, i+1) =  1;
    end
    
    % A_ineq * theta <= b_ineq
    A_ineq = D; 
    b_ineq = zeros(b, 1);

    % --- [Step 3] Optimization Setup ---
    % Initial guess
    theta0 = zeros(num_params, 1); 
    
    options = optimoptions('fmincon', ...
        'Display', 'off', ... 
        'Algorithm', 'sqp', ...
        'ConstraintTolerance', 1e-6, ...
        'StepTolerance', 1e-6, ...
        'MaxFunctionEvaluations', 50000);
        
    % Define objective function
    obj_fun = @(param) -log_likelihood_bspline(param, R_vec, Rf_vec, ...
        Basis_Precomputed, Smooth_AllR, Smooth_AllR_RND, months, use_delta, alpha, beta);

    % --- [Step 4] Run Optimization ---
    try
        [theta_hat, neg_LL, ~] = fmincon( ...
            obj_fun, theta0, A_ineq, b_ineq, [], [], [], [], [], options);
    catch
        % Fallback if heavy fail
        theta_hat = theta0;
        neg_LL = 1e10;
    end
        
    log_lik = -neg_LL;
    
    % Post-estimation call
    [~, delta_vec, M_vec] = log_likelihood_bspline(theta_hat, R_vec, Rf_vec, ...
        Basis_Precomputed, Smooth_AllR, Smooth_AllR_RND, months, use_delta, alpha, beta);
end