clear; clc;

Path_MainFolder  = 'D:\Google\我的雲端硬碟\學術｜研究與論文\論文著作\MLE Pricing Kernel';
Path_Data        = 'D:\Google\我的雲端硬碟\學術｜研究與論文\論文著作\CDI Method';
Path_Output      = fullfile(Path_MainFolder, 'Code', '11  Output');
Path_Output_Plot = fullfile(Path_MainFolder, 'Code', '11  Output - Plot');


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
Path_Code_11 = fullfile(Path_MainFolder, 'Code', ...
    '11  MLE -  Exponential Polynomial with Distortion and Monotonicity (full range)');
addpath(Path_Code_11);

Target_Points = 10002;
R_axis = generate_non_uniform_grid(Global_Min_R, Global_Max_R, Target_Points);
logR   = log(R_axis);

clear Path_Data_02 Target_TTM data input_filename year
clear col_data i this_field fields


%% Plot setting

set(groot, 'defaultAxesFontName', 'Times New Roman');
set(groot, 'defaultAxesFontSize', 12);
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultLineMarkerFaceColor','auto');


%% Plot M curve

% -------- User-specified (alpha, beta) by L --------
alpha_val = 1.00;
beta_val  = 0.90;

targets = [
    struct('L', 1, 'alpha', alpha_val, 'beta', beta_val);
    struct('L', 2, 'alpha', alpha_val, 'beta', beta_val);
    struct('L', 3, 'alpha', alpha_val, 'beta', beta_val)
];

select_rows = cell(3,1);
tol = 1e-9;
files = dir(fullfile(Path_Output, 'MLE_gamma_L_*.mat'));

for k = 1:numel(targets)
    L_i   = targets(k).L;
    a_tar = targets(k).alpha;
    b_tar = targets(k).beta;

    chosen = [];
    chosen_file = '';

    for f = 1:numel(files)
        Sfile = fullfile(Path_Output, files(f).name);
        S     = load(Sfile, 'L','alpha','beta','gamma_hat');

        if isfield(S,'L') && isfield(S,'alpha') && isfield(S,'beta')
            if S.L==L_i && abs(S.alpha-a_tar)<=tol && abs(S.beta-b_tar)<=tol
                chosen = S;
                chosen_file = files(f).name;
                break
            end
        end
    end
    if isempty(chosen)
        error('Exact file not found.');
    end
    select_rows{L_i} = table(L_i, a_tar, b_tar, string(chosen_file), ...
        'VariableNames', {'L','alpha','beta','file'});
end

% plot
figure('Position',[50 80 500 400]);
layout = tiledlayout(1, 1, 'TileSpacing', 'Compact', 'Padding', 'None');
nexttile;
hold on;

colors = get(groot, 'defaultAxesColorOrder');

for L_i = 3:3
    row = select_rows{L_i};

    Sfile = fullfile(Path_Output, row.file);
    data  = load(Sfile, 'gamma_hat', 'L', 'alpha', 'beta');
    gamma_hat = data.gamma_hat;
    L_current = data.L;
    fprintf('\nL = %d, alpha = %.2f, beta = %.2f, gamma_hat = %s\n', ...
        row.L, row.alpha, row.beta, mat2str(data.gamma_hat, 4));

    poly_sum = zeros(size(logR));
    for l = 1:L_current
        poly_sum = poly_sum + gamma_hat(l) * (logR.^l);
    end
    M_curve = exp(-poly_sum);

    plot(R_axis, M_curve, 'LineWidth', 1.8, ...
        'Color', colors(L_i,:), ...
        'DisplayName', sprintf('$L=%d$', row.L));
end

hold off;
xlabel('$R$');
ylabel('$M(R)$');
legend('show','Location','northeast','Box','off');
grid on;
xlim([0.8 1.2]);
set(gca,'LooseInset',get(gca,'TightInset'));

% Output
out_png = sprintf('plot_M_curve_alpha_%.2f_beta_%.2f.png', alpha_val, beta_val);
saveas(gcf, fullfile(Path_Output_Plot, out_png));


%% Parameter sets (six cases)

param_list = {
    struct('L',1,'alpha',0.95,'beta',0.90)
    struct('L',2,'alpha',0.95,'beta',0.90)
    struct('L',3,'alpha',1.00,'beta',0.90)
    struct('L',1,'alpha',1.00,'beta',1.00)
    struct('L',2,'alpha',1.00,'beta',1.00)
    struct('L',3,'alpha',1.00,'beta',1.00)
};


%% Compute risk preference indices and plot (given L, alpha, beta)

clc
colors = get(groot, 'defaultAxesColorOrder');

measures = {'ARA','RRA','AP','RP','AT','RT','A5','R5','A6','R6'};

% Define R_axis
Target_Points = 10002;
R_axis = generate_non_uniform_grid(Global_Min_R, Global_Max_R, Target_Points);
x      = log(R_axis);

file_cache = containers.Map();
files = dir(fullfile(Path_Output, 'MLE_gamma_L_*.mat'));

% ============================================================
%  Main Loop
% ============================================================
% for idx = 1:length(param_list)
for idx = 1:length(param_list)
    L_target     = param_list{idx}.L;
    alpha_target = param_list{idx}.alpha;
    beta_target  = param_list{idx}.beta;
    
    % ----- Find and Load file -----
    chosen_file = '';
    for f = 1:numel(files)
        Sfile = fullfile(Path_Output, files(f).name);
        S     = load(Sfile, 'L','alpha','beta','gamma_hat');
        if S.L == L_target && S.alpha == alpha_target && S.beta == beta_target
            chosen_file = files(f).name;
            gamma_hat = S.gamma_hat;
            break;
        end
    end
    if isempty(chosen_file)
        error('MLE file not found for L=%d alpha=%.2f beta=%.2f', ...
              L_target, alpha_target, beta_target);
    end
    fprintf('Processing: %s (L=%d)\n', chosen_file, L_target);
    
    % ----- Analytical Calculation -----
    % P(x) and its derivatives w.r.t x (where x = ln R)
    P              = zeros(size(x)); % P(x)
    P_prime        = zeros(size(x)); % P'(x)
    P_double_prime = zeros(size(x)); % P''(x)
    P_triple_prime = zeros(size(x)); % P'''(x)
    P_quad_prime   = zeros(size(x)); % P''''(x)
    P_penta_prime  = zeros(size(x)); % P'''''(x)
    
    for l = 1:L_target
        % P(x) = sum( gamma_l * x^l )
        P = P + gamma_hat(l) * (x.^l);
        
        % P'(x)
        P_prime = P_prime + l * gamma_hat(l) * (x.^(l-1));
        
        if l >= 2
            % P''(x)
            P_double_prime = P_double_prime + l * (l-1) * gamma_hat(l) * (x.^(l-2));
        end
        
        if l >= 3
            % P'''(x)
            P_triple_prime = P_triple_prime + l * (l-1) * (l-2) * gamma_hat(l) * (x.^(l-3));
        end
        
        if l >= 4
            % P''''(x)
            P_quad_prime = P_quad_prime + l * (l-1) * (l-2) * (l-3) * gamma_hat(l) * (x.^(l-4));
        end
        
        if l >= 5
            % P'''''(x)
            P_penta_prime = P_penta_prime + l * (l-1) * (l-2) * (l-3) * (l-4) * gamma_hat(l) * (x.^(l-5));
        end
    end
    
    % The SDF level: M(R) = exp(-P(ln R))
    % (ignoring constant exp(delta_t) as it cancels in ratios)
    M_val = exp(-P); 
    
    % --- Chain Rule Derivation for M derivatives ---
    % Let y = ln M = -P(x)
    % Derivatives of y w.r.t x:
    y_x     = -P_prime;
    y_xx    = -P_double_prime;
    y_xxx   = -P_triple_prime;
    y_xxxx  = -P_quad_prime;
    y_xxxxx = -P_penta_prime;
    
    % Derivatives of y w.r.t R (using Chain Rule):
    % y' = dy/dR = (dy/dx) * (1/R)
    y1 = y_x ./ R_axis;
    
    % y'' = d^2y/dR^2
    y2 = (y_xx - y_x) ./ (R_axis.^2);
    
    % y''' = d^3y/dR^3
    y3 = (y_xxx - 3*y_xx + 2*y_x) ./ (R_axis.^3);

    % y'''' (4th derivative of y w.r.t R)
    % Coeffs: 1, -6, 11, -6
    y4 = (y_xxxx - 6*y_xxx + 11*y_xx - 6*y_x) ./ (R_axis.^4);

    % y''''' (5th derivative of y w.r.t R)
    % Coeffs: 1, -10, 35, -50, 24
    y5 = (y_xxxxx - 10*y_xxxx + 35*y_xxx - 50*y_xx + 24*y_x) ./ (R_axis.^5);
    
    % Calculate Absolute Derivatives M1, M2, M3
    % M1 = M'
    M1 = M_val .* y1;
    
    % M2 = M''
    M2 = M_val .* (y2 + y1.^2);
    
    % M3 = M'''
    M3 = M_val .* (y3 + 3*y1.*y2 + y1.^3);

    % M4 = M'''' (For 5th Order Index)
    M4 = M_val .* (y4 + 4*y1.*y3 + 3*y2.^2 + 6*(y1.^2).*y2 + y1.^4);

    % M5 = M''''' (For 6th Order Index)
    term_5_1 = y5;
    term_5_2 = 5 * y1 .* y4;
    term_5_3 = 10 * y2 .* y3;
    term_5_4 = 10 * (y1.^2) .* y3;
    term_5_5 = 15 * y1 .* (y2.^2);
    term_5_6 = 10 * (y1.^3) .* y2;
    term_5_7 = y1.^5;
    M5 = M_val .* (term_5_1 + term_5_2 + term_5_3 + term_5_4 + term_5_5 + term_5_6 + term_5_7);
    
    % --- Risk Indices (Computed directly from M, M1, M2, M3) ---
    % 1. ARA = -M' / M
    %    RRA = R * ARA
    risk_pref.ARA = -M1 ./ M_val;    
    risk_pref.RRA = R_axis .* risk_pref.ARA;
    
    % 2. AP = -M'' / M'
    %    RP = R * AP
    risk_pref.AP  = -M2 ./ M1;    
    risk_pref.RP  = R_axis .* risk_pref.AP;
    
    % 3. AT = -M''' / M''
    %    RT = R * AT
    risk_pref.AT  = -M3 ./ M2;    
    risk_pref.RT  = R_axis .* risk_pref.AT;

    % 4. Order 5 (Edginess) -> Formula: -M'''' / M'''
    risk_pref.A5  = -M4 ./ M3;
    risk_pref.R5  = R_axis .* risk_pref.A5;
    
    % 5. Order 6 -> Formula: -M''''' / M''''
    risk_pref.A6  = -M5 ./ M4;
    risk_pref.R6  = R_axis .* risk_pref.A6;
    
    % ============================================================
    %  Output Data Table (Filtered & Formatted)
    % ============================================================
    T_full = table(R_axis, M_val, M1, M2, M3, ...
                   risk_pref.ARA, risk_pref.RRA, ...
                   risk_pref.AP, risk_pref.RP, ...
                   risk_pref.AT, risk_pref.RT, ...
                   risk_pref.A5, risk_pref.R5, ...
                   risk_pref.A6, risk_pref.R6, ...
                   'VariableNames', {'R', 'M', 'M1', 'M2', 'M3', ...
                                     'ARA', 'RRA', 'AP', 'RP', 'AT', 'RT', ...
                                     'A5', 'R5', 'A6', 'R6'});
    
    mask = (R_axis >= 1.17) & (R_axis <= 1.18);
    T_out = T_full(mask, :);
    
    varNames = T_out.Properties.VariableNames;
    for v = 1:numel(varNames)
        T_out.(varNames{v}) = round(T_out.(varNames{v}), 8);
    end
    
    alpha_str = sprintf('%.2f', alpha_target);
    beta_str  = sprintf('%.2f', beta_target);
    csv_filename = sprintf('RiskTable_L_%d_alpha_%s_beta_%s.csv', L_target, alpha_str, beta_str);
    full_csv_path = fullfile(Path_Output, csv_filename);
    
    writetable(T_out, full_csv_path);
    fprintf('Saved filtered table: %s (Rows: %d)\n', csv_filename, height(T_out));
    
    % ============================================================
    %  Plotting (Using full data, NOT filtered data)
    % ============================================================
    for m = 1:length(measures)
        key    = measures{m};
        Y_plot = risk_pref.(key);
        
        fig = figure('Position',[50 80 450 400], 'Visible', 'off');
        layout = tiledlayout(1, 1, 'TileSpacing', 'Compact', 'Padding', 'None');
        nexttile;
        hold on;
        
        plot(R_axis, Y_plot, 'LineWidth', 1.5, ...
            'Color', colors(L_target,:));
        
        grid on;
        hold off;
        xlabel('$R$', 'Interpreter', 'latex');
        ylabel(key);
        
        axis tight; 
        switch key
            case 'ARA', y_limits = [0, 4.5]; 
            case 'RRA', y_limits = [0, 3.6];
            case 'AP',  y_limits = [0, 10];
            case 'RP',  y_limits = [0, 10];
            case 'AT',  y_limits = [0, 10];
            case 'RT',  y_limits = [0, 10];
            case 'A5',  y_limits = [0, 16];
            case 'R5',  y_limits = [0, 12];
            case 'A6',  y_limits = [0, 16];
            case 'R6',  y_limits = [0, 14];
            otherwise,  axis tight; y_limits = ylim;
        end
        ylim(y_limits);
        xlim([0.8 1.2]);
        
        set(gca,'LooseInset',get(gca,'TightInset'));
        
        out_png = sprintf('%s_L_%d_alpha_%s_beta_%s.png', ...
                          key, L_target, alpha_str, beta_str);
        saveas(fig, fullfile(Path_Output_Plot, out_png));
        close(fig);
    end
end
fprintf('\nAll plots and tables completed.\n');


%% Compute delta_t and plot (given L, alpha, beta)

clc
colors = get(groot, 'defaultAxesColorOrder');
files = dir(fullfile(Path_Output, 'MLE_gamma_L_*.mat'));
use_delta = true;
dates = Realized_Return.date;
date_objs = datetime(dates, 'ConvertFrom', 'yyyymmdd');

for i = 1:length(param_list)
    L_tar = param_list{i}.L;
    a_tar = param_list{i}.alpha;
    b_tar = param_list{i}.beta;
    
    % 4.1 尋找對應的 .mat 檔案
    chosen_file = '';
    gamma_hat = [];
    tol = 1e-9;
    
    for f = 1:numel(files)
        Sfile = fullfile(Path_Output, files(f).name);
        S_info = load(Sfile, 'L', 'alpha', 'beta'); 
        if S_info.L == L_tar && abs(S_info.alpha - a_tar) < tol && abs(S_info.beta - b_tar) < tol
            chosen_file = files(f).name;
            S_full = load(Sfile, 'gamma_hat');
            gamma_hat = S_full.gamma_hat;
            break;
        end
    end
    
    if isempty(chosen_file)
        warning('Case L=%d, alpha=%.2f, beta=%.2f not found. Skipping.', L_tar, a_tar, b_tar);
        continue;
    end
    
    fprintf('Processing: L=%d, alpha=%.2f, beta=%.2f\n', L_tar, a_tar, b_tar);
    
    % 4.2 計算 Delta Time Series 
    [~, delta_vec, ~] = log_likelihood_function(gamma_hat, ...
        Realized_Return.realized_ret, Risk_Free_Rate, L_tar, ...
        Smooth_AllR, Smooth_AllR_RND, dates, use_delta, a_tar, b_tar);
    
    % 4.3 繪圖
    fig = figure('Position', [50 80 450 400], 'Visible', 'on');
    plot(date_objs, delta_vec, 'LineWidth', 1.5, 'Color', colors(L_tar,:));
    
    % X 軸設定
    datetick('x', 'yyyy', 'KeepTicks');
    xlim([min(date_objs), max(date_objs)]);
    
    % Y 軸設定
    ylabel('$\delta_t$', 'Rotation', 0);
    ylim([-0.09 0.03]);
    grid on;
    
    set(gca, 'LooseInset', get(gca, 'TightInset'));
    
    % 4.4 存檔
    out_name = sprintf('plot_delta_L_%d_alpha_%.2f_beta_%.2f.png', L_tar, a_tar, b_tar);
    saveas(fig, fullfile(Path_Output_Plot, out_name));
    % close(fig);
end

fprintf('All delta plots completed.\n');