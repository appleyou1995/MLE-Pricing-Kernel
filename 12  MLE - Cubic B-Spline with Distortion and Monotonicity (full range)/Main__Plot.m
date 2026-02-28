clear; clc;

Path_MainFolder  = 'D:\Google\我的雲端硬碟\學術｜研究與論文\論文著作\MLE Pricing Kernel';
Path_Data        = 'D:\Google\我的雲端硬碟\學術｜研究與論文\論文著作\CDI Method';
% Path_Output      = fullfile(Path_MainFolder, 'Code', '12  Output');
% Path_Output_Plot = fullfile(Path_MainFolder, 'Code', '12  Output - Plot');
Path_Output      = fullfile(Path_MainFolder, 'Code', '12  Output - Quartic');
Path_Output_Plot = fullfile(Path_MainFolder, 'Code', '12  Output - Plot (Quartic)');


%% Load the data

Target_TTM = 30;

% Load realized gross returns (R_{t+1})
Path_Data_01 = fullfile(Path_Data, 'Code', '01  輸出資料');
FileName = ['Old_Realized_Return_TTM_', num2str(Target_TTM), '.csv'];
Realized_Return = readtable(fullfile(Path_Data_01, FileName));

% Load risk-free rate R_f^t
Path_Data_01_main = fullfile(Path_Data, 'Code', '01  原始資料處理');
FileName = 'Risk_Free_Rate.csv';
Risk_Free_Rate_All = readtable(fullfile(Path_Data_01_main, FileName));
Risk_Free_Rate = Risk_Free_Rate_All.rf_gross_29d;

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


%% Parameter sets & Plot Groups

% param_list = {
%     struct('b',4,'alpha',0.95,'beta',0.90)
%     struct('b',5,'alpha',1.00,'beta',0.90)
%     struct('b',6,'alpha',1.00,'beta',0.90)
%     struct('b',7,'alpha',0.95,'beta',1.20)
%     struct('b',8,'alpha',0.95,'beta',1.20)
%     struct('b',9,'alpha',0.95,'beta',1.20)
%     struct('b',4,'alpha',1.00,'beta',1.00)
%     struct('b',5,'alpha',1.00,'beta',1.00)
%     struct('b',6,'alpha',1.00,'beta',1.00)
%     struct('b',7,'alpha',1.00,'beta',1.00)
%     struct('b',8,'alpha',1.00,'beta',1.00)
%     struct('b',9,'alpha',1.00,'beta',1.00)
% };

param_list = {
    % --- Distorted cases ---
    struct('b', 4, 'alpha', 0.95, 'beta', 0.90)
    struct('b', 5, 'alpha', 1.00, 'beta', 0.90)
    struct('b', 6, 'alpha', 1.00, 'beta', 0.90)
    struct('b', 7, 'alpha', 1.00, 'beta', 0.90)
    struct('b', 8, 'alpha', 0.95, 'beta', 1.10)
    struct('b', 9, 'alpha', 1.00, 'beta', 1.10)
    
    % --- Undistorted cases ---
    struct('b', 4, 'alpha', 1.00, 'beta', 1.00)
    struct('b', 5, 'alpha', 1.00, 'beta', 1.00)
    struct('b', 6, 'alpha', 1.00, 'beta', 1.00)
    struct('b', 7, 'alpha', 1.00, 'beta', 1.00)
    struct('b', 8, 'alpha', 1.00, 'beta', 1.00)
    struct('b', 9, 'alpha', 1.00, 'beta', 1.00)
};

% 設定繪圖群組：前六組為 Distorted，後六組為 Undistorted
plot_groups = {
    struct('indices', 1:6,  'suffix', 'distorted',   'title', 'Distorted')
    struct('indices', 7:12, 'suffix', 'undistorted', 'title', 'Undistorted')
};

tol = 1e-9;
files = dir(fullfile(Path_Output, 'MLE_BSpline_b_*.mat'));
colors = get(groot, 'defaultAxesColorOrder');


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot M curve (Distorted vs Undistorted)
for g = 1:numel(plot_groups)
    group = plot_groups{g};
    
    fig = figure('Position',[50 80 600 500]);
    set(fig, 'GraphicsSmoothing', 'on');
    set(fig, 'Renderer', 'painters');
    layout = tiledlayout(1, 1, 'TileSpacing', 'Compact', 'Padding', 'None');
    nexttile;
    hold on;
    
    for k = 1:numel(group.indices)
        idx = group.indices(k);
        p = param_list{idx};
        
        chosen_file = '';
        theta_hat = [];
        for f = 1:numel(files)
            try
                S = load(fullfile(Path_Output, files(f).name), 'b_val','alpha','beta','theta_hat');
            catch, continue; end
            
            if S.b_val == p.b && abs(S.alpha - p.alpha) <= tol && abs(S.beta - p.beta) <= tol
                chosen_file = files(f).name;
                theta_hat = S.theta_hat;
                break;
            end
        end
        
        if isempty(chosen_file), continue; end
        
        % --- Quartic B-spline 設定 ---
        n        = 4;           % Quartic
        k_order  = n + 1;    
        min_knot = Global_Min_R;
        max_knot = Global_Max_R;
        num_basis_function = p.b + 1;
        num_breaks = num_basis_function - k_order + 2;
        breaks = linspace(min_knot, max_knot, num_breaks);
        knots  = augknt(breaks, k_order);
        
        sp = spmak(knots, theta_hat');
        Spline_Val = fnval(sp, R_axis);
        M_curve = exp(Spline_Val);
        
        legend_str = sprintf('$b=%d, \\,\\alpha=%.2f, \\,\\beta=%.2f$', p.b, p.alpha, p.beta);
        plot(R_axis, M_curve, 'LineWidth', 2, ...
            'Color', colors(k,:), ... 
            'DisplayName', legend_str);
    end
    
    xlabel('$R$');
    ylabel('$M(R)$');
    grid on;
    xlim([0.8 1.2]);
    ylim([0   1.8]);
    legend('show', 'Location', 'northeast', 'Box', 'off', 'FontSize', 11);
    
    out_name = sprintf('plot_M_curve_Group_%s.png', group.suffix);
    exportgraphics(fig, fullfile(Path_Output_Plot, out_name), 'Resolution', 300);
    hold off;
end


%% Plot Delta Time Series (Distorted vs Undistorted)

dates = Realized_Return.date;
date_objs = datetime(dates, 'ConvertFrom', 'yyyymmdd');
use_delta = true;
months = Smooth_AllR.Properties.VariableNames;
T_len = length(Realized_Return.realized_ret);

for g = 1:numel(plot_groups)
    group = plot_groups{g};
    
    fig = figure('Position',[50 80 500 500]);
    set(fig, 'GraphicsSmoothing', 'on');
    set(fig, 'Renderer', 'painters');
    layout = tiledlayout(1, 1, 'TileSpacing', 'Compact', 'Padding', 'None');
    nexttile;
    hold on;
    
    for k = 1:numel(group.indices)
        idx = group.indices(k);
        p = param_list{idx};
        
        chosen_file = '';
        theta_hat = [];
        for f = 1:numel(files)
            try
                S = load(fullfile(Path_Output, files(f).name), 'b_val','alpha','beta','theta_hat');
            catch, continue; end
            
            if S.b_val == p.b && abs(S.alpha - p.alpha) <= tol && abs(S.beta - p.beta) <= tol
                chosen_file = files(f).name;
                theta_hat = S.theta_hat;
                break;
            end
        end
        
        if isempty(chosen_file), continue; end
        
        % Quartic B-spline 基礎矩陣預先計算
        n        = 4;
        k_order  = n + 1;    
        num_basis_function = p.b + 1;
        num_breaks = num_basis_function - k_order + 2;
        breaks = linspace(Global_Min_R, Global_Max_R, num_breaks);
        knots  = augknt(breaks, k_order);
        
        Basis_Precomputed = cell(T_len, 1);
        for t = 1:T_len
            col_name = months{t};
            R_axis_t = Smooth_AllR.(col_name);
            Basis_Precomputed{t} = spcol(knots, k_order, R_axis_t(:));
        end
        
        % 計算 Delta Time Series
        [~, ~, delta_vec] = log_likelihood_bspline(theta_hat, ...
            Realized_Return.realized_ret, Risk_Free_Rate, ...
            Basis_Precomputed, Smooth_AllR, Smooth_AllR_RND, ...
            months, use_delta, p.alpha, p.beta);
        
        legend_str = sprintf('$b=%d, \\,\\alpha=%.2f, \\,\\beta=%.2f$', p.b, p.alpha, p.beta);
        plot(date_objs, delta_vec, 'LineWidth', 1.5, ...
            'Color', colors(k,:), ... 
            'DisplayName', legend_str);
    end
    
    ylabel('$\delta_t$', 'Rotation', 0, 'HorizontalAlignment', 'right');
    grid on;
    xlim([min(date_objs), max(date_objs)]);
    % ylim([-0.09 0.03]);    
    xtickformat('yyyy');
    legend('show', 'Location', 'southwest', 'Box', 'off', 'FontSize', 11);
    
    out_name = sprintf('plot_delta_Group_%s.png', group.suffix);
    exportgraphics(fig, fullfile(Path_Output_Plot, out_name), 'Resolution', 300);
    hold off;
end


%% Compute and Save Risk Preference Indices (Up to Temperance)

clc
measures = {'ARA','RRA','AP','RP','AT','RT'}; 
risk_results = cell(length(param_list), 1);

for idx = 1:length(param_list)
    p = param_list{idx};
    
    chosen_file = '';
    theta_hat = [];
    for f = 1:numel(files)
        try
            S_info = load(fullfile(Path_Output, files(f).name), 'b_val','alpha','beta','theta_hat');
        catch, continue; end
        if S_info.b_val == p.b && abs(S_info.alpha - p.alpha) < tol && abs(S_info.beta - p.beta) < tol
            chosen_file = files(f).name;
            theta_hat = S_info.theta_hat;
            break;
        end
    end
    
    if isempty(chosen_file), continue; end
    
    % Quartic B-Spline Derivatives
    n        = 4; 
    k_order  = n + 1;    
    num_basis_function = p.b + 1;
    num_breaks = num_basis_function - k_order + 2;
    breaks = linspace(Global_Min_R, Global_Max_R, num_breaks);
    knots  = augknt(breaks, k_order);
    sp = spmak(knots, theta_hat');
    
    % 計算導數
    S_0 = fnval(sp, R_axis);             
    S_1 = fnval(fnder(sp, 1), R_axis);   
    S_2 = fnval(fnder(sp, 2), R_axis);   
    S_3 = fnval(fnder(sp, 3), R_axis);   
    
    M_val = exp(S_0);
    M1    = S_1 .* M_val;
    M2    = (S_2 + S_1.^2) .* M_val;
    M3    = (S_3 + 3*S_2.*S_1 + S_1.^3) .* M_val;
    
    risk_pref = struct();
    risk_pref.ARA = -S_1; 
    risk_pref.RRA = R_axis .* risk_pref.ARA;
    risk_pref.AP  = -M2 ./ M1;
    risk_pref.RP  = R_axis .* risk_pref.AP;
    risk_pref.AT  = -M3 ./ M2;
    risk_pref.RT  = R_axis .* risk_pref.AT;
    
    risk_results{idx} = risk_pref;
    
    % 輸出 Table (CSV)
    T_out = table(R_axis, M_val, M1, M2, M3, ...
                   risk_pref.ARA, risk_pref.RRA, ...
                   risk_pref.AP, risk_pref.RP, ...
                   risk_pref.AT, risk_pref.RT, ...
                   'VariableNames', {'R', 'M', 'M1', 'M2', 'M3', ...
                                     'ARA', 'RRA', 'AP', 'RP', 'AT', 'RT'});
    mask_csv = (R_axis >= 1.09) & (R_axis <= 1.11);
    T_save = T_out(mask_csv, :);
    for v = 1:numel(T_save.Properties.VariableNames)
        T_save.(T_save.Properties.VariableNames{v}) = round(T_save.(T_save.Properties.VariableNames{v}), 8);
    end
    csv_name = sprintf('RiskTable_b_%d_alpha_%.2f_beta_%.2f.csv', p.b, p.alpha, p.beta);
    writetable(T_save, fullfile(Path_Output, csv_name));
end


%% Plot Risk Preference Indices (Distorted vs Undistorted)

for m = 1:length(measures)
    measure_key = measures{m};
    
    for g = 1:numel(plot_groups)
        group = plot_groups{g};
        
        fig = figure('Position', [50 80 500 500]);
        set(fig, 'GraphicsSmoothing', 'on');
        set(fig, 'Renderer', 'painters');
        tiledlayout(1, 1, 'TileSpacing', 'Compact', 'Padding', 'None');
        nexttile;
        hold on;
        
        for k = 1:numel(group.indices)
            idx = group.indices(k);
            p = param_list{idx};
            
            if isempty(risk_results{idx}), continue; end
            
            Y_plot = risk_results{idx}.(measure_key);
            legend_str = sprintf('$b=%d, \\,\\alpha=%.2f, \\,\\beta=%.2f$', p.b, p.alpha, p.beta);
            plot(R_axis, Y_plot, 'LineWidth', 2, 'Color', colors(k,:), 'DisplayName', legend_str);
        end
        
        xlabel('$R$'); ylabel(measure_key, 'Interpreter', 'none'); grid on;
        xlim([0.8 1.2]);
        
        % 設定 Y 軸範圍 (可依據實際結果微調)
        switch measure_key
            case 'ARA', ylim([0, 5]);   case 'RRA', ylim([0, 5]);
            case 'AP',  ylim([0, 16]);  case 'RP',  ylim([0, 16]);
            case 'AT',  ylim([0, 30]);  case 'RT',  ylim([0, 30]);
        end
        
        legend('show', 'Location', 'northeast', 'Box', 'off', 'FontSize', 11);
        
        out_name = sprintf('Group_Risk_%s_%s.png', measure_key, group.suffix);
        exportgraphics(fig, fullfile(Path_Output_Plot, out_name), 'Resolution', 300);
        close(fig);
    end
end
fprintf('Risk plots & tables done.\n');


%% Compute PIT, Plot Histogram, and Save Statistics

clc
fprintf('\nGenerating PIT Histograms and Statistics...\n');
stats_list = [];
stats_cnt = 1;

for i = 1:length(param_list)
    p = param_list{i};
    chosen_file = '';
    theta_hat = [];
    
    for f = 1:numel(files)
        try
            S_info = load(fullfile(Path_Output, files(f).name), 'b_val', 'alpha', 'beta', 'theta_hat'); 
        catch, continue; end
        
        if S_info.b_val == p.b && abs(S_info.alpha - p.alpha) < tol && abs(S_info.beta - p.beta) < tol
            chosen_file = files(f).name;
            theta_hat = S_info.theta_hat;
            break;
        end
    end
    
    if isempty(chosen_file), continue; end
    
    T_len = length(Realized_Return.realized_ret);
    Basis_Precomputed = cell(T_len, 1);
    
    % Quartic B-spline 設定
    n        = 4;
    k_order  = n + 1;    
    num_basis_function = p.b + 1;
    num_breaks = num_basis_function - k_order + 2;
    breaks = linspace(Global_Min_R, Global_Max_R, num_breaks);
    knots  = augknt(breaks, k_order);
    
    for t = 1:T_len
        col_name = months{t};
        R_axis_t = Smooth_AllR.(col_name);
        Basis_Precomputed{t} = spcol(knots, k_order, R_axis_t(:));
    end
    
    % 取得 PIT 向量
    [~, ~, ~, ~, pit_vec] = log_likelihood_bspline(theta_hat, ...
        Realized_Return.realized_ret, Risk_Free_Rate, ...
        Basis_Precomputed, Smooth_AllR, Smooth_AllR_RND, ...
        months, use_delta, p.alpha, p.beta);
    
    % --- 統計檢定 ---
    pd_uniform = makedist('Uniform', 'Lower', 0, 'Upper', 1);
    [~, p_ks, ks_stat] = kstest(pit_vec, 'CDF', pd_uniform);
    
    z_norm = norminv(pit_vec); 
    z_norm = z_norm(isfinite(z_norm)); 
    [~, p_jb, jb_stat] = jbtest(z_norm);
    lags_to_test = [1, 5, 10];
    [~, p_lb, lb_stat] = lbqtest(z_norm, 'Lags', lags_to_test);
    
    stats_list(stats_cnt).b     = p.b;
    stats_list(stats_cnt).alpha = p.alpha;
    stats_list(stats_cnt).beta  = p.beta;
    stats_list(stats_cnt).KS_Stat = ks_stat;
    stats_list(stats_cnt).KS_Pval = p_ks;
    stats_list(stats_cnt).JB_Stat = jb_stat;
    stats_list(stats_cnt).JB_Pval = p_jb;
    stats_list(stats_cnt).LB_Stat_Lag1  = lb_stat(1);
    stats_list(stats_cnt).LB_Pval_Lag1  = p_lb(1);
    stats_list(stats_cnt).LB_Stat_Lag5  = lb_stat(2);
    stats_list(stats_cnt).LB_Pval_Lag5  = p_lb(2);
    stats_list(stats_cnt).LB_Stat_Lag10 = lb_stat(3);
    stats_list(stats_cnt).LB_Pval_Lag10 = p_lb(3);
    stats_cnt = stats_cnt + 1;
    
    % --- 繪製 Histogram ---
    fig = figure('Position', [50 80 450 400], 'Visible', 'on');
    layout = tiledlayout(1, 1, 'TileSpacing', 'Compact', 'Padding', 'None');
    nexttile;
    hold on;
    histogram(pit_vec, 20, 'Normalization', 'pdf', 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'w');
    yline(1, 'r--', 'LineWidth', 2, 'DisplayName', 'Uniform(0,1)');
    hold off;
    
    xlabel('$z_t$'); ylabel('Density');
    xlim([0 1]); ylim([0 1.8]); grid on;
    set(gca, 'LooseInset', get(gca, 'TightInset'));
    
    out_name = sprintf('plot_PIT_b_%d_alpha_%.2f_beta_%.2f.png', p.b, p.alpha, p.beta);
    saveas(fig, fullfile(Path_Output_Plot, out_name));
    close(fig);
end

if ~isempty(stats_list)
    T_stats = struct2table(stats_list);
    vars = T_stats.Properties.VariableNames;
    for i = 1:numel(vars)
        if isnumeric(T_stats.(vars{i}))
            T_stats.(vars{i}) = round(T_stats.(vars{i}), 4);
        end
    end
    csv_filename = fullfile(Path_Output, 'PIT_Test_Statistics_BSpline.csv');
    writetable(T_stats, csv_filename);
    fprintf('\nStatistics saved to: %s\n', csv_filename);
end
fprintf('All PIT histograms and statistics completed.\n');


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot M curve

% -------- User-specified (alpha, beta) by b (Knots) --------
alpha_val = 1.00;
beta_val  = 1.00;

b_range = 4:9;
targets = [];
for b_idx = b_range
    targets = [targets; struct('b', b_idx, 'alpha', alpha_val, 'beta', beta_val)]; %#ok<AGROW>
end

select_rows = cell(length(b_range), 1);
plot_data   = [];

tol = 1e-9;
files = dir(fullfile(Path_Output, 'MLE_BSpline_b_*.mat'));

for k = 1:numel(targets)
    b_tar     = targets(k).b;
    alpha_tar = targets(k).alpha;
    beta_tar  = targets(k).beta;
    
    chosen_file = '';
    for f = 1:numel(files)
        try
            S = load(fullfile(Path_Output, files(f).name), 'b_val','alpha','beta','theta_hat');
        catch, continue; end
        
        if S.b_val == b_tar && abs(S.alpha - alpha_tar) <= tol && abs(S.beta - beta_tar) <= tol
            chosen_file = files(f).name;
            
            % --- 計算 M Curve ---
            theta_hat = S.theta_hat;
            
            n        = 3;
            k_order  = n + 1;    
            min_knot = Global_Min_R;
            max_knot = Global_Max_R;
            num_basis_function = b_tar + 1;
            num_breaks = num_basis_function - k_order + 2;
            breaks = linspace(min_knot, max_knot, num_breaks);
            knots  = augknt(breaks, k_order);
            
            % 建立 Spline
            sp = spmak(knots, theta_hat');
            
            % 計算 Spline Value: S(R)
            Spline_Val = fnval(sp, R_axis);
            
            % M(R) = exp( S(R) )
            M_curve = exp(Spline_Val);
            
            plot_data = [plot_data; struct('b', b_tar, 'M', M_curve)];     %#ok<AGROW>
            break;
        end
    end
    if isempty(chosen_file)
        warning('File not found for b=%d', b_tar);
    end
end

% Plot
figure('Position',[50 80 500 400]);
layout = tiledlayout(1, 1, 'TileSpacing', 'Compact', 'Padding', 'None');
nexttile;
hold on;
colors = get(groot, 'defaultAxesColorOrder');

for k = 1:length(plot_data)
    D = plot_data(k);
    plot(R_axis, D.M, 'LineWidth', 1.8, ...
        'Color', colors(k,:), ...
        'DisplayName', sprintf('$b=%d$', D.b));
end

hold off;
xlabel('$R$');
ylabel('$M(R)$');
legend('show','Location','northeast','Box','off');
grid on;
xlim([0.8 1.2]);
set(gca,'LooseInset',get(gca,'TightInset'));

% Output
out_png = sprintf('plot_M_curve_BSpline_alpha_%.2f_beta_%.2f.png', alpha_val, beta_val);
saveas(gcf, fullfile(Path_Output_Plot, out_png));


%% Compute risk preference indices and plot

clc
colors = get(groot, 'defaultAxesColorOrder');
measures = {'ARA','RRA','AP','RP','AT','RT','A5','R5','A6','R6'}; 

files = dir(fullfile(Path_Output, 'MLE_BSpline_b_*.mat'));

% Main Loop
for idx = 1:length(param_list)
    b_target     = param_list{idx}.b;
    alpha_target = param_list{idx}.alpha;
    beta_target  = param_list{idx}.beta;
    
    % ----- Find and Load -----
    chosen_file = '';
    theta_hat = [];
    
    for f = 1:numel(files)
        try
            S_info = load(fullfile(Path_Output, files(f).name), 'b_val','alpha','beta');
        catch, continue; end
        
        if S_info.b_val == b_target && ...
           abs(S_info.alpha - alpha_target) < tol && ...
           abs(S_info.beta - beta_target) < tol
            
            chosen_file = files(f).name;
            S_full = load(fullfile(Path_Output, files(f).name), 'theta_hat');
            theta_hat = S_full.theta_hat;
            break;
        end
    end
    
    if isempty(chosen_file)
        fprintf('File not found for b=%d, alpha=%.2f, beta=%.2f. Skipping...\n', ...
            b_target, alpha_target, beta_target);
        continue;
    end
    fprintf('Processing: %s (b=%d)\n', chosen_file, b_target);
    
    % ----- B-Spline Derivatives -----
    % Reconstruct Knots
    n        = 3;
    k_order  = n + 1;    
    min_knot = Global_Min_R;
    max_knot = Global_Max_R;
    num_basis_function = b_target + 1;
    num_breaks = num_basis_function - k_order + 2;
    breaks = linspace(min_knot, max_knot, num_breaks);
    knots  = augknt(breaks, k_order);
    
    % Create Spline Structure
    sp = spmak(knots, theta_hat');
    
    % Compute Spline Values and Derivatives w.r.t R 
    S_0 = fnval(sp, R_axis);             % S(R)
    S_1 = fnval(fnder(sp, 1), R_axis);   % S'(R)
    S_2 = fnval(fnder(sp, 2), R_axis);   % S''(R)
    S_3 = fnval(fnder(sp, 3), R_axis);   % S'''(R)
    
    % Calculate M Derivatives
    % Model: M(R) = exp( S(R) )  (Ignoring delta)
    % M'   = S' * M
    % M''  = (S'' + (S')^2) * M
    % M''' = (S''' + 3*S''*S' + (S')^3) * M
    
    M_val = exp(S_0);
    M1    = S_1 .* M_val;
    M2    = (S_2 + S_1.^2) .* M_val;
    M3    = (S_3 + 3*S_2.*S_1 + S_1.^3) .* M_val;

    % --- M4 (假設 S4=0) ---
    Bracket_4 = 4.*S_3.*S_1 + 3.*S_2.^2 + 6.*S_2.*S_1.^2 + S_1.^4;
    M4        = Bracket_4 .* M_val;

    % --- M5 (假設 S4=0, S5=0) ---
    Bracket_5 = 10.*S_3.*S_2 + 10.*S_3.*S_1.^2 + 15.*S_1.*S_2.^2 + ...
                10.*(S_1.^3).*S_2 + S_1.^5;
    M5        = Bracket_5 .* M_val;
    
    % --- Risk Indices ---
    % ARA = -M' / M = -S'
    risk_pref.ARA = -S_1; 
    risk_pref.RRA = R_axis .* risk_pref.ARA;
    
    % AP = -M'' / M' = -(S'' + (S')^2) / S' = -S''/S' - S'
    risk_pref.AP  = -M2 ./ M1;
    risk_pref.RP  = R_axis .* risk_pref.AP;
    
    % AT = -M''' / M'' 
    risk_pref.AT  = -M3 ./ M2;
    risk_pref.RT  = R_axis .* risk_pref.AT;

    % Order 5: -M4 / M3
    risk_pref.A5  = -M4 ./ M3;
    risk_pref.R5  = R_axis .* risk_pref.A5;
    
    % Order 6: -M5 / M4
    risk_pref.A6  = -M5 ./ M4;
    risk_pref.R6    = R_axis .* risk_pref.A6;
    
    % Output Table
    T_out = table(R_axis, M_val, M1, M2, M3, ...
                   risk_pref.ARA, risk_pref.RRA, ...
                   risk_pref.AP, risk_pref.RP, ...
                   risk_pref.AT, risk_pref.RT, ...
                   risk_pref.A5, risk_pref.R5, ...
                   risk_pref.A6, risk_pref.R6, ...
                   'VariableNames', {'R', 'M', 'M1', 'M2', 'M3', ...
                                     'ARA', 'RRA', 'AP', 'RP', 'AT', 'RT', ...
                                     'A5', 'R5', 'A6', 'R6'});
    
    % Filter for saving CSV
    mask_csv = (R_axis >= 1.09) & (R_axis <= 1.11);
    T_save = T_out(mask_csv, :);

    varNames = T_save.Properties.VariableNames;
    for v = 1:numel(varNames)
        T_save.(varNames{v}) = round(T_save.(varNames{v}), 8);
    end
    
    csv_name = sprintf('RiskTable_b_%d_alpha_%.2f_beta_%.2f.csv', b_target, alpha_target, beta_target);
    writetable(T_save, fullfile(Path_Output, csv_name));

    fprintf('Saved risk table: %s\n', csv_name);
    
    % ============================================================
    %  Trend Check Table (0.8, 1.0, 1.2)
    % ============================================================
    target_R_points = [0.8, 1.0, 1.2];
    trend_data = struct();
    trend_data.Measure = measures(:); 
    
    for k = 1:length(target_R_points)
        r_val = target_R_points(k);
        raw_name = sprintf('Val_at_%.1f', r_val);
        col_name = strrep(raw_name, '.', '_');
        
        vals = zeros(length(measures), 1);
        for m = 1:length(measures)
            key = measures{m};
            y_data = risk_pref.(key);
            vals(m) = interp1(R_axis, y_data, r_val, 'pchip');
        end
        trend_data.(col_name) = vals;
    end
    
    T_trend = struct2table(trend_data);
    
    % Format
    varNamesTrend = T_trend.Properties.VariableNames;
    for v = 1:numel(varNamesTrend)
        if isnumeric(T_trend.(varNamesTrend{v}))
            T_trend.(varNamesTrend{v}) = round(T_trend.(varNamesTrend{v}), 4);
        end
    end
    
    trend_csv_name = sprintf('TrendCheck_b_%d_alpha_%.2f_beta_%.2f.csv', b_target, alpha_target, beta_target);
    writetable(T_trend, fullfile(Path_Output, trend_csv_name));
    fprintf('Saved trend check table: %s\n', trend_csv_name);
    
    % Plotting
    for m = 1:length(measures)
        key = measures{m};
        Y_plot = risk_pref.(key);
        
        fig = figure('Position',[50 80 450 400], 'Visible', 'off');
        layout = tiledlayout(1, 1, 'TileSpacing', 'Compact', 'Padding', 'None');
        nexttile;
        hold on;
        
        plot(R_axis, Y_plot, 'LineWidth', 1.5, 'Color', colors(b_target-3,:));
        
        grid on; hold off;
        xlabel('$R$', 'Interpreter', 'latex');
        ylabel(key);
        axis tight;
        switch key
            case 'ARA', y_limits = [0, 5]; 
            case 'RRA', y_limits = [0, 5];
            case 'AP',  y_limits = [0, 16];
            case 'RP',  y_limits = [0, 16];
            case 'AT',  y_limits = [0, 30];
            case 'RT',  y_limits = [0, 30];
            case 'A5',  y_limits = [0, 12];
            case 'R5',  y_limits = [0, 10];
            case 'A6',  y_limits = [0, 16];
            case 'R6',  y_limits = [0, 14];
            otherwise,  axis tight; y_limits = ylim;
        end
        ylim(y_limits);
        xlim([0.8 1.2]);
        
        out_png = sprintf('%s_b_%d_alpha_%.2f_beta_%.2f.png', key, b_target, alpha_target, beta_target);
        saveas(fig, fullfile(Path_Output_Plot, out_png));
        close(fig);
    end
end

fprintf('Risk plots & tables done.\n');


%% Compute delta_t and plot (given b, alpha, beta)

clc
colors = get(groot, 'defaultAxesColorOrder');
files = dir(fullfile(Path_Output, 'MLE_BSpline_b_*.mat'));
use_delta = true;

months = Smooth_AllR.Properties.VariableNames;
dates  = Realized_Return.date;
date_objs = datetime(dates, 'ConvertFrom', 'yyyymmdd');
tol = 1e-9;

for i = 1:length(param_list)
    b_tar     = param_list{i}.b;
    alpha_tar = param_list{i}.alpha;
    beta_tar  = param_list{i}.beta;
    
    chosen_file = '';
    theta_hat = [];
    tol = 1e-9;
    
    for f = 1:numel(files)
        Sfile = fullfile(Path_Output, files(f).name);
        try
            S_info = load(Sfile, 'b_val', 'alpha', 'beta'); 
        catch
            continue;
        end
        
        if S_info.b_val == b_tar && ...
           abs(S_info.alpha - alpha_tar) < tol && ...
           abs(S_info.beta - beta_tar) < tol
       
            chosen_file = files(f).name;
            S_full = load(Sfile, 'theta_hat');
            theta_hat = S_full.theta_hat;
            break;
        end
    end
    
    if isempty(chosen_file)
        warning('Case b=%d, alpha=%.2f, beta=%.2f not found. Skipping.', b_tar, alpha_tar, beta_tar);
        continue;
    end
    
    fprintf('Processing: b=%d, alpha=%.2f, beta=%.2f\n', b_tar, alpha_tar, beta_tar);
    
    T = length(Realized_Return.realized_ret);
    Basis_Precomputed = cell(T, 1);
    
    % Knots
    n        = 3;
    k_order  = n + 1;    
    min_knot = Global_Min_R;
    max_knot = Global_Max_R;
    num_basis_function = b_tar + 1;
    num_breaks = num_basis_function - k_order + 2;
    breaks = linspace(min_knot, max_knot, num_breaks);
    knots  = augknt(breaks, k_order);
    
    % Basis Matrix
    for t = 1:T
        col_name = months{t};
        R_axis_t = Smooth_AllR.(col_name);
        R_axis_t = R_axis_t(:);
        Basis_Precomputed{t} = spcol(knots, k_order, R_axis_t);
    end

    % --- Delta Time Series ---
    [~, ~, delta_vec] = log_likelihood_bspline(theta_hat, ...
        Realized_Return.realized_ret, Risk_Free_Rate, ...
        Basis_Precomputed, ...
        Smooth_AllR, Smooth_AllR_RND, ...
        months, ...
        use_delta, alpha_tar, beta_tar);
    
    % --- Plot ---
    fig = figure('Position', [50 80 450 400], 'Visible', 'on');
    
    color_idx = b_tar - 3;
    plot(date_objs, delta_vec, 'LineWidth', 1.5, 'Color', colors(color_idx,:));
    
    datetick('x', 'yyyy', 'KeepTicks');
    xlim([min(date_objs), max(date_objs)]);
    
    ylabel('$\delta_t$', 'Rotation', 0);
    % ylim([0 1.56]);
    grid on;
    
    set(gca, 'LooseInset', get(gca, 'TightInset'));
    
    out_name = sprintf('plot_delta_b_%d_alpha_%.2f_beta_%.2f.png', b_tar, alpha_tar, beta_tar);
    saveas(fig, fullfile(Path_Output_Plot, out_name));
    % close(fig);
end

fprintf('All delta plots completed.\n');


%% Compute PIT, Plot Histogram, and Save Statistics

clc
fprintf('\nGenerating PIT Histograms and Statistics...\n');
colors = get(groot, 'defaultAxesColorOrder');
files = dir(fullfile(Path_Output, 'MLE_BSpline_b_*.mat'));
use_delta = true;
months = Smooth_AllR.Properties.VariableNames;
dates = Realized_Return.date;

stats_list = [];
stats_cnt = 1;

for i = 1:length(param_list)
    b_tar     = param_list{i}.b;
    alpha_tar = param_list{i}.alpha;
    beta_tar  = param_list{i}.beta;
    
    chosen_file = '';
    theta_hat = [];
    
    for f = 1:numel(files)
        Sfile = fullfile(Path_Output, files(f).name);
        try
            S_info = load(Sfile, 'b_val', 'alpha', 'beta'); 
        catch, continue; end
        
        if S_info.b_val == b_tar && ...
           abs(S_info.alpha - alpha_tar) < tol && ...
           abs(S_info.beta - beta_tar) < tol
       
            chosen_file = files(f).name;
            S_full = load(Sfile, 'theta_hat');
            theta_hat = S_full.theta_hat;
            break;
        end
    end
    
    if isempty(chosen_file)
        warning('Case b=%d, alpha=%.2f, beta=%.2f not found. Skipping.', b_tar, alpha_tar, beta_tar);
        continue;
    end
    
    fprintf('Processing PIT: b=%d, alpha=%.2f, beta=%.2f\n', b_tar, alpha_tar, beta_tar);
    
    % Basis Matrix Calculation (Need to re-compute for PIT)
    T = length(Realized_Return.realized_ret);
    Basis_Precomputed = cell(T, 1);
    
    n        = 3;
    k_order  = n + 1;    
    min_knot = Global_Min_R;
    max_knot = Global_Max_R;
    num_basis_function = b_tar + 1;
    num_breaks = num_basis_function - k_order + 2;
    breaks = linspace(min_knot, max_knot, num_breaks);
    knots  = augknt(breaks, k_order);
    
    for t = 1:T
        col_name = months{t};
        R_axis_t = Smooth_AllR.(col_name);
        R_axis_t = R_axis_t(:);
        Basis_Precomputed{t} = spcol(knots, k_order, R_axis_t);
    end

    % 計算 PIT 向量
    [~, ~, ~, ~, pit_vec] = log_likelihood_bspline(theta_hat, ...
        Realized_Return.realized_ret, Risk_Free_Rate, ...
        Basis_Precomputed, ...
        Smooth_AllR, Smooth_AllR_RND, ...
        months, ...
        use_delta, alpha_tar, beta_tar);
    
    % ============================================================
    %  Statistical Tests
    % ============================================================
    
    % KS Test
    pd_uniform = makedist('Uniform', 'Lower', 0, 'Upper', 1);
    [~, p_ks, ks_stat] = kstest(pit_vec, 'CDF', pd_uniform);
    
    % Berkowitz
    z_norm = norminv(pit_vec); 
    z_norm = z_norm(isfinite(z_norm)); 
    
    % JB Test
    [~, p_jb, jb_stat] = jbtest(z_norm);
    
    % LB Test
    lags_to_test = [1, 5, 10];
    [~, p_lb, lb_stat] = lbqtest(z_norm, 'Lags', lags_to_test);
    
    % --- 儲存結果 ---
    stats_list(stats_cnt).b     = b_tar;
    stats_list(stats_cnt).alpha = alpha_tar;
    stats_list(stats_cnt).beta  = beta_tar;
    
    stats_list(stats_cnt).KS_Stat = ks_stat;
    stats_list(stats_cnt).KS_Pval = p_ks;
    stats_list(stats_cnt).JB_Stat = jb_stat;
    stats_list(stats_cnt).JB_Pval = p_jb;
    
    stats_list(stats_cnt).LB_Stat_Lag1  = lb_stat(1);
    stats_list(stats_cnt).LB_Pval_Lag1  = p_lb(1);
    stats_list(stats_cnt).LB_Stat_Lag5  = lb_stat(2);
    stats_list(stats_cnt).LB_Pval_Lag5  = p_lb(2);
    stats_list(stats_cnt).LB_Stat_Lag10 = lb_stat(3);
    stats_list(stats_cnt).LB_Pval_Lag10 = p_lb(3);
    
    stats_cnt = stats_cnt + 1;

    % ============================================================
    %  Plotting Histogram
    % ============================================================
    fig = figure('Position', [50 80 450 400], 'Visible', 'on');
    layout = tiledlayout(1, 1, 'TileSpacing', 'Compact', 'Padding', 'None');
    nexttile;
    hold on;
    
    h = histogram(pit_vec, 20, 'Normalization', 'pdf', ...
        'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'w');
    
    yline(1, 'r--', 'LineWidth', 2, 'DisplayName', 'Uniform(0,1)');
    
    hold off;
    
    xlabel('$z_t$');
    ylabel('Density');
    xlim([0 1]);
    ylim([0 1.8]);
    
    grid on;
    set(gca, 'LooseInset', get(gca, 'TightInset'));
    
    out_name = sprintf('plot_PIT_b_%d_alpha_%.2f_beta_%.2f.png', b_tar, alpha_tar, beta_tar);
    saveas(fig, fullfile(Path_Output_Plot, out_name));
    % close(fig);
end

% Save Statistics
if ~isempty(stats_list)
    T_stats = struct2table(stats_list);
    vars = T_stats.Properties.VariableNames;
    for i = 1:numel(vars)
        if isnumeric(T_stats.(vars{i}))
            T_stats.(vars{i}) = round(T_stats.(vars{i}), 4);
        end
    end
    csv_filename = fullfile(Path_Output, 'PIT_Test_Statistics_BSpline.csv');
    writetable(T_stats, csv_filename);
    fprintf('\nStatistics saved to: %s\n', csv_filename);
end

fprintf('All PIT histograms and statistics completed.\n');