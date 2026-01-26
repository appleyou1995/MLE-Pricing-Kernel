clear; clc;

Path_MainFolder  = 'D:\Google\我的雲端硬碟\學術｜研究與論文\論文著作\MLE Pricing Kernel';
Path_Data        = 'D:\Google\我的雲端硬碟\學術｜研究與論文\論文著作\CDI Method';
Path_Output      = fullfile(Path_MainFolder, 'Code', '12  Output');
Path_Output_Plot = fullfile(Path_MainFolder, 'Code', '12  Output - Plot');


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


%% Plot M curve

% -------- User-specified (alpha, beta) by b (Knots) --------
alpha_val = 1.00;
beta_val  = 1.00;

targets = [
    struct('b', 4, 'alpha', alpha_val, 'beta', beta_val);
    struct('b', 6, 'alpha', alpha_val, 'beta', beta_val);
    struct('b', 8, 'alpha', alpha_val, 'beta', beta_val)
];

select_rows = cell(3,1);
plot_data   = [];

tol = 1e-9;
files = dir(fullfile(Path_Output, 'MLE_BSpline_b_*.mat'));

for k = 1:numel(targets)
    beta_tar = targets(k).b;
    alpha_tar = targets(k).alpha;
    beta_tar = targets(k).beta;
    
    chosen_file = '';
    for f = 1:numel(files)
        try
            S = load(fullfile(Path_Output, files(f).name), 'b_val','alpha','beta','theta_hat');
        catch, continue; end
        
        if S.b_val == beta_tar && abs(S.alpha - alpha_tar) <= tol && abs(S.beta - beta_tar) <= tol
            chosen_file = files(f).name;
            
            % --- 計算 M Curve ---
            theta_hat = S.theta_hat;
            
            % 重建 Knots (邏輯需與估計程式一致)
            n         = 3; 
            num_knots = n + beta_tar + 2;
            knots = linspace(Global_Min_R, Global_Max_R, num_knots);
            knots(1:(n+1))     = Global_Min_R;
            knots((end-n):end) = Global_Max_R;
            
            % 建立 Spline
            sp = spmak(knots, theta_hat');
            
            % 計算 Spline Value: S(R)
            Spline_Val = fnval(sp, R_axis);
            
            % M(R) = exp( S(R) )  <-- 根據 Eq(4)，指數是正的
            % 注意：這裡忽略 delta (常數項)，只畫形狀
            M_curve = exp(Spline_Val);
            
            % 存起來畫圖用
            plot_data = [plot_data; struct('b', beta_tar, 'M', M_curve)]; %#ok<AGROW>
            break;
        end
    end
    if isempty(chosen_file)
        warning('File not found for b=%d', beta_tar);
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


%% Parameter sets

param_list = {
    struct('b',4,'alpha',0.95,'beta',0.90)
    struct('b',6,'alpha',1.00,'beta',0.90)
    struct('b',8,'alpha',1.00,'beta',0.90)
    struct('b',4,'alpha',1.00,'beta',1.00)
    struct('b',6,'alpha',1.00,'beta',1.00)
    struct('b',8,'alpha',1.00,'beta',1.00)
};


%% Compute risk preference indices and plot

clc
colors = get(groot, 'defaultAxesColorOrder');

% 建議移除 A5, R5, A6, R6 (Cubic Spline 4階以上微分為0)
measures = {'ARA','RRA','AP','RP','AT','RT','A5','R5','A6','R6'}; 

files = dir(fullfile(Path_Output, 'MLE_BSpline_b_*.mat'));

% Main Loop
for idx = 1:length(param_list)
    b_target     = param_list{idx}.b;
    alpha_target = param_list{idx}.alpha;
    beta_target  = param_list{idx}.beta;
    
    % ----- Find and Load -----
    chosen_file = '';
    for f = 1:numel(files)
        try
            S = load(fullfile(Path_Output, files(f).name), 'b_val','alpha','beta','theta_hat');
        catch, continue; end
        
        if S.b_val == b_target && abs(S.alpha - alpha_target)<1e-9 && abs(S.beta - beta_target)<1e-9
            chosen_file = files(f).name;
            theta_hat = S.theta_hat;
            break;
        end
    end
    
    if isempty(chosen_file)
        fprintf('File not found for b=%d, skipping...\n', b_target);
        continue;
    end
    fprintf('Processing: %s (b=%d)\n', chosen_file, b_target);
    
    % ----- B-Spline Derivatives -----
    % Reconstruct Knots
    n         = 3;
    num_knots = n + b_target + 2;
    knots = linspace(Global_Min_R, Global_Max_R, num_knots);
    knots(1:(n+1))     = Global_Min_R;
    knots((end-n):end) = Global_Max_R;
    
    % Create Spline Structure
    sp = spmak(knots, theta_hat');
    
    % Compute Spline Values and Derivatives w.r.t R
    % S(R), S'(R), S''(R), S'''(R)
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
    risk_pref.A5 = -M4 ./ M3;
    risk_pref.R5 = R_axis .* risk_pref.A5;
    
    % Order 6: -M5 / M4
    risk_pref.A6 = -M5 ./ M4;
    risk_pref.R6 = R_axis .* risk_pref.A6;
    
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
    
    % Plotting
    for m = 1:length(measures)
        key = measures{m};
        Y_plot = risk_pref.(key);
        
        fig = figure('Position',[50 80 450 400], 'Visible', 'off');
        layout = tiledlayout(1, 1, 'TileSpacing', 'Compact', 'Padding', 'None');
        nexttile;
        hold on;
        
        plot(R_axis, Y_plot, 'LineWidth', 1.5, 'Color', colors(b_target/2-1,:));
        
        grid on; hold off;
        xlabel('$R$', 'Interpreter', 'latex');
        ylabel(key);
        axis tight;
        switch key
            case 'ARA', y_limits = [0, 5]; 
            case 'RRA', y_limits = [0, 5];
            case 'AP',  y_limits = [0, 16];
            case 'RP',  y_limits = [0, 16];
            case 'AT',  y_limits = [0, 10];
            case 'RT',  y_limits = [0, 10];
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
fprintf('Done.\n');


%% Compute delta_t and plot (given b, alpha, beta)

clc
colors = get(groot, 'defaultAxesColorOrder');
files = dir(fullfile(Path_Output, 'MLE_BSpline_b_*.mat'));
use_delta = true;

months = Smooth_AllR.Properties.VariableNames;
dates  = Realized_Return.date;
date_objs = datetime(dates, 'ConvertFrom', 'yyyymmdd');

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
        
        if S_info.b_val == b_tar && abs(S_info.alpha - alpha_tar) < tol && abs(S_info.beta - beta_tar) < tol
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
    n         = 3;
    k_order   = n + 1;
    num_knots = n + b_tar + 2;
    min_knot  = Global_Min_R;
    max_knot  = Global_Max_R;
    
    knots = linspace(min_knot, max_knot, num_knots);
    knots(1:(n+1))      = min_knot; 
    knots((end-n):end)  = max_knot;
    
    % Basis Matrix
    for t = 1:T
        col_name = months{t};
        R_axis_t = Smooth_AllR.(col_name);
        R_axis_t = R_axis_t(:);
        Basis_Precomputed{t} = spcol(knots, k_order, R_axis_t);
    end

    % --- Delta Time Series ---
    [~, delta_vec, ~] = log_likelihood_bspline(theta_hat, ...
        Realized_Return.realized_ret, Risk_Free_Rate, ...
        Basis_Precomputed, ...
        Smooth_AllR, Smooth_AllR_RND, ...
        months, ...
        use_delta, alpha_tar, beta_tar);
    
    % --- Plot ---
    fig = figure('Position', [50 80 450 400], 'Visible', 'on');
    
    color_idx = b_tar/2 - 1;
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