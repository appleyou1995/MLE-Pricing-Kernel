clear; clc;

Path_MainFolder  = 'D:\Google\我的雲端硬碟\學術｜研究與論文\論文著作\MLE Pricing Kernel';
Path_Data        = 'D:\Google\我的雲端硬碟\學術｜研究與論文\論文著作\CDI Method';
Path_Output      = fullfile(Path_MainFolder, 'Code', '15  Output - without and wide');
Path_Output_plot = fullfile(Path_MainFolder, 'Code', '15  Output - without and wide (alpha groups)');


%% Load the data

Target_TTM = 30;

% Load realized gross returns (R_{t+1})
Path_Data_01 = fullfile(Path_Data, 'Code', '01  輸出資料');
FileName = ['Old_Realized_Return_TTM_', num2str(Target_TTM), '.csv'];
Realized_Return = readtable(fullfile(Path_Data_01, FileName));
R_vec = Realized_Return.realized_ret(:);

% Load risk-free rate R_f^t
Path_Data_01_main = fullfile(Path_Data, 'Code', '01  原始資料處理');
FileName = 'Risk_Free_Rate.csv';
Risk_Free_Rate_All = readtable(fullfile(Path_Data_01_main, FileName));
Risk_Free_Rate = Risk_Free_Rate_All.rf_gross_29d;
Rf_vec = Risk_Free_Rate(:);

% Load Q-measure PDF tables: R axis and corresponding f^*_t(R)
Path_Data_02 = fullfile(Path_Data, 'Code', '02  輸出資料');
Smooth_AllR     = [];
Smooth_AllR_RND = [];

years_to_merge = 1996:2021;
for year = years_to_merge
    input_filename = fullfile(Path_Data_02, sprintf('TTM_%d_RND_Tables_%d.mat', Target_TTM, year));
    if exist(input_filename, 'file')
        data = load(input_filename);
        Smooth_AllR     = [Smooth_AllR, data.Table_Smooth_AllR];           %#ok<AGROW>
        Smooth_AllR_RND = [Smooth_AllR_RND, data.Table_Smooth_AllR_RND];   %#ok<AGROW>
    else
        warning('File %s does not exist.', input_filename);
    end
end

% Find the range of R_axis
Global_Min_R = 100; 
Global_Max_R = 0;
fields = Smooth_AllR.Properties.VariableNames; 

for i = 1:numel(fields)
    this_field = fields{i};
    col_data = Smooth_AllR.(this_field); 
    
    if isnumeric(col_data) && ~isempty(col_data)
        Global_Min_R = min(Global_Min_R, min(col_data));
        Global_Max_R = max(Global_Max_R, max(col_data));
    end
end
Global_Min_R = Global_Min_R * 0.9; 
Global_Max_R = Global_Max_R * 1.1;

% Define R_axis
Target_Points = 10000;
R_axis = linspace(Global_Min_R, Global_Max_R, Target_Points)';

clear Path_Data_02 Target_TTM data input_filename year
clear col_data i this_field fields


%% Plot setting

set(groot, 'defaultAxesFontName', 'Times New Roman');
set(groot, 'defaultAxesFontSize', 12);
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultLineMarkerFaceColor','auto');


%% Parameter & Plot settings

alpha_range = 0.7:0.1:1.2;
beta_range  = 0.7:0.1:1.2;

% Quintic
n_degree = 5;
b_val = 6;
k_order = n_degree + 1;
tol = 1e-9;
files = dir(fullfile(Path_Output, 'MLE_BSpline_b_*.mat'));

% 顏色漸層設定：
num_betas = numel(beta_range);
custom_colors = turbo(num_betas);


%% Plot M curve (以 Alpha 分組，每張圖疊加所有 Beta)

% Add paths
Path_Code_15 = fullfile(Path_MainFolder, 'Code', ...
    '15  MLE - B-Spline (Monotonically decreasing constraint)');
addpath(Path_Code_15);

% === 預先計算 Basis_Original ===
fprintf('正在預計算原始資料的 B-spline 基底... ');
num_basis_function = b_val + 1;
num_breaks = num_basis_function - k_order + 2;
breaks = linspace(Global_Min_R, Global_Max_R, num_breaks);
knots  = augknt(breaks, k_order);

months = Smooth_AllR.Properties.VariableNames;
T_len = length(Realized_Return.realized_ret);
Basis_Original = cell(T_len, 1);
for t = 1:T_len
    col_name = months{t};
    R_axis_t = Smooth_AllR.(col_name);
    Basis_Original{t} = spcol(knots, k_order, R_axis_t(:));
end
B_plot = spcol(knots, k_order, R_axis);
fprintf('完成。\n');

R_vec_numeric = Realized_Return.realized_ret;
Rf_vec_numeric = Risk_Free_Rate;

% ==================== 主繪圖迴圈 ====================
% 外層迴圈：每一個 alpha 畫一張獨立的圖
for a_idx = 1:numel(alpha_range)
    current_alpha = alpha_range(a_idx);
    
    % 建立畫布
    fig = figure('Position',[50 80 600 500], 'Visible', 'on');
    set(fig, 'GraphicsSmoothing', 'on');
    set(fig, 'Renderer', 'painters');
    layout = tiledlayout(1, 1, 'TileSpacing', 'Compact', 'Padding', 'None');
    nexttile;
    hold on;
    
    plot_counter = 0;
    
    % 內層迴圈：在目前的 alpha 底下，把所有的 beta 曲線疊上來
    for b_idx = 1:numel(beta_range)
        current_beta = beta_range(b_idx);
        
        % --- 1. 尋找對應的 .mat 檔案 ---
        chosen_file = '';
        theta_hat = [];
        for f = 1:numel(files)
            try
                S = load(fullfile(Path_Output, files(f).name), 'b_val','alpha','beta','theta_hat');
            catch, continue; end
            
            if S.b_val == b_val && abs(S.alpha - current_alpha) <= tol && abs(S.beta - current_beta) <= tol
                chosen_file = files(f).name;
                theta_hat = S.theta_hat(:);
                break;
            end
        end
        
        % 如果這個參數組合在資料夾中找不到（沒跑出結果），就跳過
        if isempty(chosen_file), continue; end
        plot_counter = plot_counter + 1;
        
        % --- 2. 獲取該組合每一期的 delta_vec ---
        [~, ~, delta_vec] = log_likelihood_bspline(theta_hat, ...
            R_vec_numeric, Rf_vec_numeric, ...
            Basis_Original, Smooth_AllR, Smooth_AllR_RND, ...
            months, true, current_alpha, current_beta);
        
        % --- 3. 計算並映射到高解析度 R_axis ---
        Spline_Sum_Plot = B_plot * theta_hat; 
        Spline_Sum_Plot = max(min(Spline_Sum_Plot, 60), -60);
        
        % 矩陣廣播計算時間序列平均 SDF
        M_all = exp(delta_vec + Spline_Sum_Plot'); 
        M_curve_avg = mean(M_all, 1);
        
        % --- 4. 繪製曲線 ---
        legend_str = sprintf('$\\beta = %.2f$', current_beta);
        plot(R_axis, M_curve_avg, 'LineWidth', 1.8, ...
            'Color', custom_colors(b_idx, :), ... 
            'DisplayName', legend_str);
    end
    
    % --- 5. 圖表裝飾與優化 ---
    if plot_counter > 0
        xlabel('$R_{t+1}$');
        ylabel('$\mathrm{E}[M_t(R_{t+1})]$', 'Interpreter', 'latex');
        title(sprintf('Fixed $\\alpha = %.2f$', current_alpha), 'Interpreter', 'latex');
        grid on;
        xlim([0.8 1.2]);
        legend('show', 'Location', 'northeast', 'Box', 'off', 'FontSize', 10);
        
        % 檔案命名：
        out_name = sprintf('sdf_sensitivity_unweighted_alpha_%.2f.png', current_alpha);
        exportgraphics(fig, fullfile(Path_Output_plot, out_name), 'Resolution', 300);
        fprintf('匯出圖表: %s (共包含 %d 條 Beta 曲線)\n', out_name, plot_counter);
    end
    
    close(fig);
end
fprintf('==== 所有圖繪製完成！ ====\n');





%% Plot m curve (以 Alpha 分組，每張圖疊加所有 Beta)

% Add paths
Path_Code_15 = fullfile(Path_MainFolder, 'Code', ...
    '15  MLE - B-Spline (Monotonically decreasing constraint)');
addpath(Path_Code_15);

% === 預先計算 Basis_Original ===
fprintf('正在預計算原始資料的 B-spline 基底... ');
num_basis_function = b_val + 1;
num_breaks = num_basis_function - k_order + 2;
breaks = linspace(Global_Min_R, Global_Max_R, num_breaks);
knots  = augknt(breaks, k_order);

months = Smooth_AllR.Properties.VariableNames;
T_len = length(Realized_Return.realized_ret);
Basis_Original = cell(T_len, 1);
for t = 1:T_len
    col_name = months{t};
    R_axis_t = Smooth_AllR.(col_name);
    Basis_Original{t} = spcol(knots, k_order, R_axis_t(:));
end
B_plot = spcol(knots, k_order, R_axis);
fprintf('完成。\n');

R_vec_numeric = Realized_Return.realized_ret;
Rf_vec_numeric = Risk_Free_Rate;

% ==================== 主繪圖迴圈 ====================
% 外層迴圈：每一個 alpha 畫一張獨立的圖
for a_idx = 1:numel(alpha_range)
    current_alpha = alpha_range(a_idx);
    
    % 建立畫布
    fig = figure('Position',[50 80 600 500], 'Visible', 'on');
    set(fig, 'GraphicsSmoothing', 'on');
    set(fig, 'Renderer', 'painters');
    layout = tiledlayout(1, 1, 'TileSpacing', 'Compact', 'Padding', 'None');
    nexttile;
    hold on;
    
    plot_counter = 0;
    
    % 內層迴圈：在目前的 alpha 底下，把所有的 beta 曲線疊上來
    for b_idx = 1:numel(beta_range)
        current_beta = beta_range(b_idx);
        
        % --- 1. 尋找對應的 .mat 檔案 ---
        chosen_file = '';
        theta_hat = [];
        for f = 1:numel(files)
            try
                S = load(fullfile(Path_Output, files(f).name), 'b_val','alpha','beta','theta_hat');
            catch, continue; end
            
            if S.b_val == b_val && abs(S.alpha - current_alpha) <= tol && abs(S.beta - current_beta) <= tol
                chosen_file = files(f).name;
                theta_hat = S.theta_hat(:);
                break;
            end
        end
        
        % 如果這個參數組合在資料夾中找不到（沒跑出結果），就跳過
        if isempty(chosen_file), continue; end
        plot_counter = plot_counter + 1;
        
        % --- 2. 獲取該組合每一期的 delta_vec ---
        [~, ~, delta_vec] = log_likelihood_bspline(theta_hat, ...
            R_vec_numeric, Rf_vec_numeric, ...
            Basis_Original, Smooth_AllR, Smooth_AllR_RND, ...
            months, true, current_alpha, current_beta);
        
        % --- 3. 計算並映射到高解析度 R_axis ---
        Spline_Sum_Plot = B_plot * theta_hat; 
        Spline_Sum_Plot = max(min(Spline_Sum_Plot, 60), -60);
        
        % 計算大 M (純邊際效用)
        M_all = exp(delta_vec + Spline_Sum_Plot'); % size: [T_len, Target_Points]
        
        % 準備存放小 m (包含機率扭曲) 的矩陣
        m_all = zeros(T_len, length(R_axis));
        
        % 逐期計算機率扭曲 D'(F^P)
        for t = 1:T_len
            % 抓取第 t 期的 R_axis 與 Q測度機率密度 (RND)
            col_name = months{t};
            R_axis_t = Smooth_AllR.(col_name);
            RND_t    = Smooth_AllR_RND.(col_name);
            
            % 將 RND 內插到高解析度、統整的 R_axis 上 (超出範圍補 0)
            RND_interp = interp1(R_axis_t, RND_t, R_axis, 'linear', 0);
            
            % 取出第 t 期的 M_t 與 Rf
            M_t  = M_all(t, :)';
            Rf_t = Rf_vec_numeric(t);
            
            % 步驟 A：計算積分 I_t(x) = \int_0^x (R_f^{-1} / M_t) * f^Q dx
            integrand = (1 / Rf_t) .* (RND_interp(:) ./ M_t);
            I_t = cumtrapz(R_axis, integrand); 
            
            % 確保機率值在合理範圍 (避免後續 log(0) 或負數產生複數)
            I_t = max(min(I_t, 1 - 1e-8), 1e-8);
            
            % 步驟 B：透過反函數推導 F^P
            % F^P = exp( -(1/beta) * (-log(I_t))^(1/alpha) )
            F_P = exp( -(1/current_beta) * (-log(I_t)).^(1/current_alpha) );
            F_P = max(min(F_P, 1 - 1e-8), 1e-8);
            
            % 步驟 C：計算機率扭曲函數的一階導數 D'(F^P)
            % D'(x) = exp(-(-beta*log(x))^alpha) * (alpha*beta/x) * (-beta*log(x))^(alpha-1)
            % 這裡的 exp(...) 剛好就是 I_t
            term1 = I_t;
            term2 = (current_alpha * current_beta) ./ F_P;
            term3 = (-current_beta * log(F_P)).^(current_alpha - 1);
            D_prime = term1 .* term2 .* term3;
            
            % 步驟 D：組合得到小 m
            m_all(t, :) = M_t .* D_prime;
        end
        
        % 計算時間序列平均的小 m 曲線
        m_curve_avg = mean(m_all, 1);
        
        % --- 4. 繪製曲線 ---
        legend_str = sprintf('$\\beta = %.2f$', current_beta);
        plot(R_axis, m_curve_avg, 'LineWidth', 1.8, ...
            'Color', custom_colors(b_idx, :), ... 
            'DisplayName', legend_str);
    end
    
    % --- 5. 圖表裝飾與優化 ---
    if plot_counter > 0
        xlabel('$R_{t+1}$');
        ylabel('$\mathrm{E}[m_t(R_{t+1})]$', 'Interpreter', 'latex');
        title(sprintf('Fixed $\\alpha = %.2f$', current_alpha), 'Interpreter', 'latex');
        grid on;
        xlim([0.8 1.2]);
        legend('show', 'Location', 'northeast', 'Box', 'off', 'FontSize', 10);
        
        % 檔案命名：
        out_name = sprintf('sdf_sensitivity_weighted_alpha_%.2f.png', current_alpha);
        exportgraphics(fig, fullfile(Path_Output_plot, out_name), 'Resolution', 300);
        fprintf('匯出圖表: %s (共包含 %d 條 Beta 曲線)\n', out_name, plot_counter);
    end
    
    close(fig);
end
fprintf('==== 所有圖繪製完成！ ====\n');