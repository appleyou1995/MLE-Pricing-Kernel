clear; clc; close all;

Path_MainFolder  = 'D:\Google\我的雲端硬碟\學術｜研究與論文\論文著作\MLE Pricing Kernel';
Path_Data        = 'D:\Google\我的雲端硬碟\學術｜研究與論文\論文著作\CDI Method';

Path_Code        = fullfile(Path_MainFolder, 'Code');
Path_Code_16     = fullfile(Path_MainFolder, 'Code', '16  MLE - B-Spline (Extension)');


%% Setting

set(groot, 'defaultAxesFontName', 'Times New Roman');
set(groot, 'defaultAxesFontSize', 12);
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');

R_axis = linspace(0.7, 1.3, 5000)';
b_val = 6;


%% Plot Groups

% --- Group 1: Benchmark vs Non-Recession ---
Group1.Name = 'G1_Recession';
Group1.OutPath = fullfile(Path_Code, '16  Output (Non-Recession Periods) - Plot');
Group1.YLim_Set = struct('ARA', [0, 5],  'RRA', [0, 5], ...
                         'AP',  [4, 12], 'RP',  [4, 12], ...
                         'AT',  [4, 12], 'RT',  [4, 12]);
Group1.Lines = {
    struct('Name', 'Benchmark',    'TTM', 30, 'Mode', 'Full',          'a_u', 1.00, 'b_u', 1.00, 'a_d', 1.00, 'b_d', 0.90, 'LS', '--')
    struct('Name', 'Non-Recession','TTM', 30, 'Mode', 'Non_Recession', 'a_u', 1.00, 'b_u', 1.00, 'a_d', 0.95, 'b_d', 0.95, 'LS', '-')
};

% --- Group 2: Benchmark vs Volatility ---
Group2.Name = 'G2_Volatility';
Group2.OutPath = fullfile(Path_Code, '16  Output (High vs. Low Volatility) - Plot');
Group2.YLim_Set = struct('ARA', [0, 12],  'RRA', [0, 12], ...
                         'AP',  [-5, 25], 'RP',  [-5, 25], ...
                         'AT',  [-5, 25], 'RT',  [-5, 25]);
Group2.Lines = {
    struct('Name', 'Benchmark',       'TTM', 30, 'Mode', 'Full',     'a_u', 1.00, 'b_u', 1.00, 'a_d', 1.00, 'b_d', 0.90, 'LS', '--')
    struct('Name', 'High Volatility', 'TTM', 30, 'Mode', 'High_Vol', 'a_u', 1.00, 'b_u', 1.00, 'a_d', 0.90, 'b_d', 1.10, 'LS', '-')
    struct('Name', 'Low Volatility',  'TTM', 30, 'Mode', 'Low_Vol',  'a_u', 1.00, 'b_u', 1.00, 'a_d', 0.90, 'b_d', 1.15, 'LS', '-')
};

% --- Group 3: Benchmark vs Term Structure (TTM) ---
Group3.Name = 'G3_TermStructure';
Group3.OutPath = fullfile(Path_Code, '16  Output (TTM Comparison) - Plot');
Group3.YLim_Set = struct('ARA', [0, 5],  'RRA', [0, 5], ...
                         'AP',  [0, 12], 'RP',  [0, 12], ...
                         'AT',  [-2, 12], 'RT',  [-2, 12]);
Group3.Lines = {
    struct('Name', 'TTM=30 (BM)', 'TTM', 30,  'Mode', 'Full', 'a_u', 1.00, 'b_u', 1.00, 'a_d', 1.00, 'b_d', 0.90, 'LS', '--')
    struct('Name', 'TTM=60',      'TTM', 60,  'Mode', 'Full', 'a_u', 1.00, 'b_u', 1.00, 'a_d', 0.95, 'b_d', 0.85, 'LS', '-')
    struct('Name', 'TTM=90',      'TTM', 90,  'Mode', 'Full', 'a_u', 1.00, 'b_u', 1.00, 'a_d', 1.00, 'b_d', 0.85, 'LS', '-')
    struct('Name', 'TTM=180',     'TTM', 180, 'Mode', 'Full', 'a_u', 1.00, 'b_u', 1.00, 'a_d', 0.90, 'b_d', 0.95, 'LS', '-')
};


%% Group selecting

% All_Groups = {Group1, Group2, Group3};
All_Groups = {Group3};
Status_List = {'Undistorted', 'Distorted'};
Risk_Metrics = {'ARA','RRA','AP','RP','AT','RT'};
colors = lines(7);


%% Main Plotting Loop

clc; close all;

for g = 1:numel(All_Groups)
    grp = All_Groups{g};
    fprintf('\nProcessing Group: %s\n', grp.Name);
    
    for s = 1:numel(Status_List)
        status = Status_List{s};
        is_distorted = strcmp(status, 'Distorted');
        
        % 初始化 M(R) 圖片
        fig_M = figure('Name', sprintf('%s_M_%s', grp.Name, status), 'Position', [50 50 500 400]);
        tiledlayout(fig_M, 1, 1, 'TileSpacing', 'Compact', 'Padding', 'None');
        nexttile;
        hold on;
        
        % 初始化 Risk Indices 圖片
        figs_R = gobjects(numel(Risk_Metrics), 1);
        for rm = 1:numel(Risk_Metrics)
            figs_R(rm) = figure('Name', sprintf('%s_%s_%s', grp.Name, Risk_Metrics{rm}, status), 'Position', [50 50 500 400]); 
            set(figs_R(rm), 'GraphicsSmoothing', 'on');
            set(figs_R(rm), 'Renderer', 'painters');
            tiledlayout(figs_R(rm), 1, 1, 'TileSpacing', 'Compact', 'Padding', 'None');
            nexttile;
            hold on;
        end
        
        % 開始畫這個 Group 的每一條線
        for L = 1:numel(grp.Lines)
            line_info = grp.Lines{L};
            
            % 根據目前的狀態決定取 a_d/b_d 還是 a_u/b_u
            alpha_val = sum([line_info.a_d, line_info.a_u] .* [is_distorted, ~is_distorted]);
            beta_val  = sum([line_info.b_d, line_info.b_u] .* [is_distorted, ~is_distorted]);
            
            % 1. 取得 Theta Hat (自動跨資料夾搜尋)
            theta_hat = get_theta(Path_Code, line_info.Mode, line_info.TTM, alpha_val, beta_val);
            if isempty(theta_hat), continue; end
            
            % 2. 載入這條線專屬的 Knots 並計算 M 與 Risk
            bounds_file = fullfile(Path_Code_16, sprintf('Global_Bounds_TTM_%d.mat', line_info.TTM));
            [M_curve, Risk_Struct] = calc_spline_metrics(theta_hat, bounds_file, b_val, R_axis);
            
            % Legend 文字
            leg_str = sprintf('%s ($\\alpha=%.2f, \\beta=%.2f$)', line_info.Name, alpha_val, beta_val);
            
            % --- 畫 M(R) ---
            figure(fig_M);
            plot(R_axis, M_curve, line_info.LS, 'LineWidth', 2, 'Color', colors(L,:), 'DisplayName', leg_str);
            
            % --- 畫 Risk Indices ---
            for rm = 1:numel(Risk_Metrics)
                figure(figs_R(rm));
                plot(R_axis, Risk_Struct.(Risk_Metrics{rm}), line_info.LS, 'LineWidth', 2, 'Color', colors(L,:), 'DisplayName', leg_str);
            end
        end
        
        % --- 圖片收尾與存檔 (M curve) ---
        figure(fig_M);
        grid on;
        xlabel('$R_{t+1}$');
        ylabel('$M_{t}(R_{t+1})$');
        xlim([0.8 1.2]);
        legend('show', 'Location', 'northeast', 'Box', 'off', 'FontSize', 10);
        exportgraphics(fig_M, fullfile(grp.OutPath, sprintf('%s_M_%s.png', grp.Name, status)), 'Resolution', 300);
        close(fig_M);
        
        % --- 圖片收尾與存檔 (Risk Indices) ---
        for rm = 1:numel(Risk_Metrics)
            figure(figs_R(rm));
            grid on;
            xlabel('$R_{t+1}$');
            ylabel(Risk_Metrics{rm}, 'Interpreter', 'none');
            xlim([0.8 1.2]);
            
            % Y 軸範圍
            current_metric = Risk_Metrics{rm};
            if isfield(grp.YLim_Set, current_metric)
                ylim(grp.YLim_Set.(current_metric));
            end
            
            legend('show', 'Location', 'northeast', 'Box', 'off', 'FontSize', 10);
            exportgraphics(figs_R(rm), ...
                fullfile(grp.OutPath, sprintf('%s_%s_%s.png', grp.Name, Risk_Metrics{rm}, status)), ...
                'Resolution', 300);
            close(figs_R(rm));
        end
    end
end
disp('All plots generated successfully!');


%% Local Functions

function theta = get_theta(Path_Code, Mode, TTM, alpha, beta)
    % 依照 Extension_Mode 的寫法，明確定義各個 Mode 的資料夾路徑
    switch Mode
        case 'Full'
            folder_path = fullfile(Path_Code, '16  Output (TTM Comparison)');
            
        case 'Non_Recession'
            folder_path = fullfile(Path_Code, '16  Output (Non-Recession Periods)');
            
        case 'High_Vol'
            folder_path = fullfile(Path_Code, '16  Output (High vs. Low Volatility)');
            
        case 'Low_Vol'
            folder_path = fullfile(Path_Code, '16  Output (High vs. Low Volatility)');
            
        otherwise
            error('未知的 Mode: %s', Mode);
    end
    
    % 定義檔名並組合完整路徑
    pattern = sprintf('MLE_BSpline_TTM_%d_%s_alpha_%.2f_beta_%.2f.mat', TTM, Mode, alpha, beta);
    file_path = fullfile(folder_path, pattern);
    
    % 讀取檔案
    if exist(file_path, 'file')
        data = load(file_path, 'theta_hat');
        theta = data.theta_hat;
    else
        warning('找不到檔案: %s\n路徑: %s', pattern, file_path);
        theta = [];
    end
end


function [M_curve, Risk] = calc_spline_metrics(theta_hat, bounds_file, b_val, R_axis)
    % 載入這條線專屬的 Global Bounds
    B = load(bounds_file, 'Global_Min_R', 'Global_Max_R');
    
    n_degree = 5;
    k_order = n_degree + 1;
    num_breaks = (b_val + 1) - k_order + 2;
    breaks = linspace(B.Global_Min_R, B.Global_Max_R, num_breaks);
    knots  = augknt(breaks, k_order);
    
    sp = spmak(knots, theta_hat');
    
    S_0 = fnval(sp, R_axis);             
    S_1 = fnval(fnder(sp, 1), R_axis);   
    S_2 = fnval(fnder(sp, 2), R_axis);   
    S_3 = fnval(fnder(sp, 3), R_axis);   
    
    M_curve = exp(S_0);
    M1 = S_1 .* M_curve;
    M2 = (S_2 + S_1.^2) .* M_curve;
    M3 = (S_3 + 3*S_2.*S_1 + S_1.^3) .* M_curve;
    
    Risk.ARA = -S_1; 
    Risk.RRA = R_axis .* Risk.ARA;
    Risk.AP  = -M2 ./ M1;
    Risk.RP  = R_axis .* Risk.AP;
    Risk.AT  = -M3 ./ M2;
    Risk.RT  = R_axis .* Risk.AT;
end