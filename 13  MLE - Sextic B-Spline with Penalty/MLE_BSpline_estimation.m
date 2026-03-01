%% Main Function: MLE Theta Estimation for Penalized B-Spline Pricing Kernel

function [theta_hat, log_lik, BIC, delta_vec, M_vec, pit_vec] = MLE_BSpline_estimation( ...
    Smooth_AllR, Smooth_AllR_RND, Realized_Return, Risk_Free_Rate, ...
    b, use_delta, alpha, beta, Global_Min_R, Global_Max_R, lambda, theta_init)

    % Settings
    rng(0);
    R_vec  = Realized_Return.realized_ret;
    Rf_vec = Risk_Free_Rate;
    months = Smooth_AllR.Properties.VariableNames;
    T      = length(R_vec);

    % --- [Step 1] Construct Knots & Precompute Basis ---
    n_degree = 6;
    k_order  = n_degree + 1;    
    min_knot = Global_Min_R;
    max_knot = Global_Max_R;

    num_basis_function = b + 1;
    num_breaks = num_basis_function - k_order + 2;
    
    if num_breaks < 2
        error('基底函數個數 b=%d 不足於支撐 degree n=%d 的 B-spline 模型。需增加 b 的值。', b, n_degree);
    end

    breaks = linspace(min_knot, max_knot, num_breaks);
    knots  = augknt(breaks, k_order);
    
    % Precompute B-Spline Basis for each day to speed up fmincon
    Basis_Precomputed = cell(T, 1);
    
    for t = 1:T
        try
            col_name = months{t};
            R_axis = Smooth_AllR.(col_name);
            R_axis = R_axis(:);
            
            B = spcol(knots, k_order, R_axis);            
            Basis_Precomputed{t} = B;
            
        catch ME
            warning('Error in precompute at t=%d (%s): %s', t, months{t}, ME.message);
            Basis_Precomputed{t} = [];
        end
    end
    
    % --- [Step 1.5] Construct Penalty Matrix (P-spline) ---
    % ★ 根據 Schneider et al. 建構懲罰矩陣 P = eta * P3 + P4
    num_params = b + 1;
    if num_params < 5
        error('基底函數數量太少，無法計算四階差分。需設定較大的 b (例如 b=20)。');
    end
    
    I_mat = eye(num_params);
    D3 = diff(I_mat, 3);  % 三階差分矩陣
    D4 = diff(I_mat, 4);  % 四階差分矩陣
    
    % 設定 eta 使得 P4 主導平滑過程，P3 僅用於抑制邊界碎波
    eta = 0.2; 
    Penalty_Matrix = eta * (D3' * D3) + (D4' * D4);

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
    if nargin > 11 && ~isempty(theta_init)
        theta0 = theta_init(:); % 使用傳入的熱啟動值
    else
        theta0 = zeros(num_params, 1); % 冷啟動
    end
    
    options = optimoptions('fmincon', ...
        'Display', 'off', ... 
        'Algorithm', 'interior-point', ...
        'OptimalityTolerance', 1e-4, ...
        'ConstraintTolerance', 1e-4, ...
        'StepTolerance', 1e-4, ...
        'MaxFunctionEvaluations', 10000);
        
    % ★ Define objective function (改用下方定義的包裝函數，加入 Penalty 與 lambda)
    obj_fun = @(param) objective_function_penalized(param, R_vec, Rf_vec, ...
        Basis_Precomputed, Smooth_AllR, Smooth_AllR_RND, months, use_delta, alpha, beta, ...
        Penalty_Matrix, lambda);
        
    % --- [Step 4] Run Optimization ---
    try
        [theta_hat, ~, ~] = fmincon( ...
            obj_fun, theta0, A_ineq, b_ineq, [], [], [], [], [], options);
    catch
        % Fallback if heavy fail
        theta_hat = theta0;
    end
        
    % ★ Post-estimation call: 使用算出的 theta_hat 計算真實的 (無懲罰) Log-likelihood
    if nargout > 3
        [log_lik_raw, BIC, delta_vec, M_vec, pit_vec] = log_likelihood_bspline(theta_hat, R_vec, Rf_vec, ...
            Basis_Precomputed, Smooth_AllR, Smooth_AllR_RND, months, use_delta, alpha, beta);
    else
        [log_lik_raw, BIC] = log_likelihood_bspline(theta_hat, R_vec, Rf_vec, ...
            Basis_Precomputed, Smooth_AllR, Smooth_AllR_RND, months, use_delta, alpha, beta);
        delta_vec = [];  M_vec = [];  pit_vec = [];
    end
    
    % 回傳真實的對數概似值
    log_lik = log_lik_raw;
end


%% Objective Function Wrapper for Penalized MLE
%  這是一個輔助函數，用來計算 fmincon 需要最小化的目標值 (-LL + lambda * Penalty)

function Neg_Penalized_LL = objective_function_penalized(theta, R_vec, Rf_vec, ...
    Basis_Precomputed, Smooth_AllR, Smooth_AllR_RND, months, use_delta, alpha, beta, ...
    Penalty_Matrix, lambda)
    
    % 1. 計算原本的純粹 Log-Likelihood
    [LL, ~, ~, ~, ~] = log_likelihood_bspline(theta, R_vec, Rf_vec, ...
        Basis_Precomputed, Smooth_AllR, Smooth_AllR_RND, months, use_delta, alpha, beta);
    
    % 2. 計算 Penalty 值: theta' * P * theta
    theta = theta(:);
    Penalty_Val = theta' * Penalty_Matrix * theta;
    
    % 3. 計算帶有懲罰項的目標函數 (因為 fmincon 是找最小值，所以 LL 要加負號)
    Neg_Penalized_LL = -LL + lambda * Penalty_Val;
end