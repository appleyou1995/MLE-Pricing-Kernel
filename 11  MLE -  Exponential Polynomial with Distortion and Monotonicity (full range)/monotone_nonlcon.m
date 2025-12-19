function [c, ceq] = monotone_nonlcon(gamma, L, R_min, R_max)
    
    % 1. 建立檢查點 (Grid Points)
    % 直接在全域範圍切 100 個點
    num_check_points = 100;
    R_grid = linspace(R_min, R_max, num_check_points)'; 
    x_grid = log(R_grid); % 轉換成 log-return 空間
    
    % 2. 計算多項式的「斜率」 (Derivative of Polynomial)
    % 向量化寫法：
    % P(x) = sum( gamma_l * x^l )
    % P'(x) = sum( l * gamma_l * x^(l-1) )
    
    % 產生次數向量 [1, 2, ..., L]
    powers = 1:L; 
    
    % 核心運算：矩陣相乘
    % Basis Matrix: 每一行是 l * x^(l-1)
    % 這裡利用 broadcasting: x_grid .^ (powers - 1)
    Basis = x_grid .^ (powers - 1); 
    
    % 對 Basis 的每一行乘上對應的係數 l (微分產生的係數)
    Basis_Derivative = Basis .* powers;
    
    % 計算每個 grid 點上的斜率 (Slope)
    % Slope 是一個 (num_check_points x 1) 的向量
    Slope = Basis_Derivative * gamma; 
    
    % 3. 設定限制式 (Inequality Constraints)
    % fmincon 要求 c(x) <= 0
    % 我們要求 Slope >= 0 (即 Slope >= tol)
    % 轉換為： tol - Slope <= 0
    
    tol = 1e-6; 
    c = tol - Slope; % 如果 Slope < tol，c 就會變成正的 -> 違反限制
    ceq = [];
end