clear; clc;

Path_MainFolder = 'D:\Google\我的雲端硬碟\學術｜研究與論文\論文著作\MLE Pricing Kernel';
Path_Data       = fullfile(Path_MainFolder, 'Code', '00  Output');
Path_Output     = fullfile(Path_MainFolder, 'Code', '17  Output - fine');


%% Load the data

Target_TTM = 30;

% Load realized gross returns (R_{t+1})
FileName = ['Realized_Return_TTM_', num2str(Target_TTM), '.csv'];
Realized_Return = readtable(fullfile(Path_Data, FileName));

% Load risk-free rate R_f^t
FileName = 'Risk_Free_GrossFactor_ByTargetTTM.csv';
Risk_Free_Rate_All = readtable(fullfile(Path_Data, FileName));

% Load Q-measure PDF tables: R axis and corresponding f^*_t(R)
Path_RND = fullfile(Path_MainFolder, 'Code', '01  Output');
Smooth_AllR = [];
Smooth_AllR_RND = [];

years_to_merge = 1996:2025;
for year = years_to_merge
    input_filename = fullfile(Path_RND, sprintf('TTM_%d_RND_Tables_%d.mat', Target_TTM, year));
    if exist(input_filename, 'file')
        data = load(input_filename);
        Smooth_AllR = [Smooth_AllR, data.Table_Smooth_AllR];               %#ok<AGROW>
        Smooth_AllR_RND = [Smooth_AllR_RND, data.Table_Smooth_AllR_RND];   %#ok<AGROW>
    else
        warning('File %s does not exist.', input_filename);
    end
end

clear Path_RND
clear data FileName input_filename year


%% Align RND, Realized Return, and Risk-Free Rate by date

% 1. RND dates from table variable names
rnd_var_names = string(Smooth_AllR_RND.Properties.VariableNames);
rnd_dates = str2double(regexp(rnd_var_names, '\d+', 'match', 'once'));

% 2. Realized return dates
realized_dates = double(Realized_Return.date);

% 3. Risk-free rate dates and values
rf_date_col = ['date_', num2str(Target_TTM)];
rf_rate_col = ['rf_gross_TTM', num2str(Target_TTM)];
rf_dates_all  = double(Risk_Free_Rate_All.(rf_date_col));
rf_values_all = double(Risk_Free_Rate_All.(rf_rate_col));
idx_rf_valid = isfinite(rf_dates_all) & isfinite(rf_values_all);
rf_dates_all  = rf_dates_all(idx_rf_valid);
rf_values_all = rf_values_all(idx_rf_valid);

% 4. Use realized return dates as the master sample
idx_realized_has_rnd = ismember(realized_dates, rnd_dates);
idx_realized_has_rf  = ismember(realized_dates, rf_dates_all);
idx_master = idx_realized_has_rnd & idx_realized_has_rf;
common_dates = realized_dates(idx_master);
common_dates = sort(common_dates);
if isempty(common_dates)
    error('No common dates across RND, Realized_Return, and Risk_Free_Rate.');
end

% 5. Locate corresponding RND columns, realized rows, and RF rows
[tf_rnd, idx_rnd] = ismember(common_dates, rnd_dates);
[tf_ret, idx_ret] = ismember(common_dates, realized_dates);
[tf_rf,  idx_rf]  = ismember(common_dates, rf_dates_all);

if ~all(tf_rnd) || ~all(tf_ret) || ~all(tf_rf)
    error('Date alignment failed.');
end

fprintf('\nDate alignment summary:\n');
fprintf('Number of RND months            : %d\n', numel(rnd_dates));
fprintf('Number of realized-return months: %d\n', numel(realized_dates));
fprintf('Number of RF months             : %d\n', numel(rf_dates_all));
fprintf('Number of aligned months        : %d\n', numel(common_dates));
fprintf('First aligned date              : %08.0f\n', common_dates(1));
fprintf('Last aligned date               : %08.0f\n', common_dates(end));

% 6. Subset and reorder all data by common dates
Smooth_AllR     = Smooth_AllR(:, idx_rnd);
Smooth_AllR_RND = Smooth_AllR_RND(:, idx_rnd);
Realized_Return = Realized_Return(idx_ret, :);
Risk_Free_Rate  = rf_values_all(idx_rf);
R_vec           = Realized_Return.realized_ret(:);
Rf_vec          = Risk_Free_Rate(:);

clear idx_realized_has_rnd idx_realized_has_rf idx_master idx_rf_valid
clear idx_rnd idx_ret idx_rf tf_rnd tf_ret tf_rf
clear common_dates rnd_dates rnd_var_names
clear rf_values_all rf_dates_all rf_rate_col rf_date_col Risk_Free_Rate_All


%% Find the range of R_axis

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

% Quintic
n_degree = 5;
param_list = {
    % --- Distorted cases ---
    struct('b', 6, 'alpha', 1.00, 'beta', 0.92)

    % --- Undistorted cases ---
    struct('b', 6, 'alpha', 1.00, 'beta', 1.00)
};



% 設定繪圖群組：前半為 Distorted，後半組為 Undistorted
plot_groups = {
    struct('indices', 1, 'suffix', 'distorted',   'title', 'Distorted')
    struct('indices', 2, 'suffix', 'undistorted', 'title', 'Undistorted')
};

k_order = n_degree + 1;
tol = 1e-9;
files = dir(fullfile(Path_Output, 'MLE_BSpline_b_*.mat'));
colors = get(groot, 'defaultAxesColorOrder');


%% Plot M curve (Distorted vs Undistorted)

% Add paths
Path_Code_15 = fullfile(Path_MainFolder, 'Code', ...
    '15  MLE - B-Spline (Monotonically decreasing constraint)');
addpath(Path_Code_15);

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
        
        % --- 1. 讀取對應的 theta_hat ---
        chosen_file = '';
        theta_hat = [];
        for f = 1:numel(files)
            try
                S = load(fullfile(Path_Output, files(f).name), 'b_val','alpha','beta','theta_hat');
            catch, continue; end
            
            if S.b_val == p.b && abs(S.alpha - p.alpha) <= tol && abs(S.beta - p.beta) <= tol
                chosen_file = files(f).name;
                theta_hat = S.theta_hat(:);
                break;
            end
        end
        if isempty(chosen_file), continue; end
        
        R_vec_numeric = Realized_Return.realized_ret;
        Rf_vec_numeric = Risk_Free_Rate;

        % --- 2. B-spline 節點設定 ---  
        num_basis_function = p.b + 1;
        num_breaks = num_basis_function - k_order + 2;
        breaks = linspace(Global_Min_R, Global_Max_R, num_breaks);
        knots  = augknt(breaks, k_order);
        
        % --- 3. 準備原始網格的 Basis (為了精確計算每一期的 delta) ---
        months = Smooth_AllR.Properties.VariableNames;
        T_len = length(Realized_Return.realized_ret);
        Basis_Original = cell(T_len, 1);
        for t = 1:T_len
            col_name = months{t};
            R_axis_t = Smooth_AllR.(col_name);
            Basis_Original{t} = spcol(knots, k_order, R_axis_t(:));
        end
        
        % --- 4. 獲取每一期的 delta_vec ---
        [~, ~, delta_vec] = log_likelihood_bspline(theta_hat, ...
            R_vec_numeric, Rf_vec_numeric, ...
            Basis_Original, Smooth_AllR, Smooth_AllR_RND, ...
            months, true, p.alpha, p.beta);
        
        % --- 5. 在統一的高解析度 R_axis 上計算 Q(R) 並裝回 delta ---
        B_plot = spcol(knots, k_order, R_axis);
        Spline_Sum_Plot = B_plot * theta_hat; % 這就是 Q(R)
        
        % 數值溢位保護 (與核心函數一致)
        Spline_Sum_Plot = max(min(Spline_Sum_Plot, 60), -60);
        
        % 矩陣廣播運算：
        % delta_vec 是 [T x 1], Spline_Sum_Plot' 是 [1 x 10000]
        % M_all 得到 [T x 10000]
        M_all = exp(delta_vec + Spline_Sum_Plot'); 
        
        % --- 6. 取時間平均並繪圖 ---
        M_curve_avg = mean(M_all, 1);
        
        legend_str = sprintf('$\\alpha=%.2f, \\,\\beta=%.2f$', p.alpha, p.beta);
        plot(R_axis, M_curve_avg, 'LineWidth', 2, ...
            'Color', colors(k,:), ... 
            'DisplayName', legend_str);
    end
    
    xlabel('$R_{t+1}$');
    % ylabel('$E[M_t(R_{t+1})]$');
    ylabel('$\mathrm{E}[M_t(R_{t+1})]$', 'Interpreter', 'latex');
    grid on;
    xlim([0.8 1.2]);
    legend('show', 'Location', 'northeast', 'Box', 'off', 'FontSize', 11);
    
    out_name = sprintf('plot_M_MarginalUtility_curve_Group_%s.png', group.suffix);
    exportgraphics(fig, fullfile(Path_Output, out_name), 'Resolution', 300);
    hold off;
end


%% Plot small-m curve: m_t(R) = M_t(R) * D'(F^P(R))
%
% Polkovnichenko and Zhao (2013):
%
%   m_t(R) = u'(R) Z(P)
%
% In our notation:
%
%   m_t(R) = M_t(R) D'(F^P(R))
%
% This section is appended after the original M-curve plotting code.
% The original code does not need to be modified.

% -------------------------------------------------------------------------
% 1. Estimate the physical CDF F^P(R) using realized returns
% -------------------------------------------------------------------------

R_realized_for_cdf = Realized_Return.realized_ret(:);
R_realized_for_cdf = R_realized_for_cdf(isfinite(R_realized_for_cdf));

if isempty(R_realized_for_cdf)
    error('No finite realized returns are available for estimating F^P(R).');
end

% Smooth kernel estimate of the physical CDF
F_P_plot = ksdensity( ...
    R_realized_for_cdf, ...
    R_axis, ...
    'Function', 'cdf');

F_P_plot = F_P_plot(:);

% Numerical protection:
% D'(x) contains log(x) and division by x.
cdf_epsilon = 1e-8;
F_P_plot = min(max(F_P_plot, cdf_epsilon), 1 - cdf_epsilon);


% -------------------------------------------------------------------------
% 2. Plot the small-m curve for each parameter group
% -------------------------------------------------------------------------

for g = 1:numel(plot_groups)

    group = plot_groups{g};

    fig = figure('Position', [50 80 600 500]);
    set(fig, 'GraphicsSmoothing', 'on');
    set(fig, 'Renderer', 'painters');

    layout = tiledlayout(1, 1, ...
        'TileSpacing', 'Compact', ...
        'Padding', 'None');

    nexttile;
    hold on;

    for k = 1:numel(group.indices)

        idx = group.indices(k);
        p = param_list{idx};

        % -----------------------------------------------------------------
        % 2.1 Load the corresponding theta_hat
        % -----------------------------------------------------------------

        chosen_file = '';
        theta_hat = [];

        for f = 1:numel(files)

            try
                S = load( ...
                    fullfile(Path_Output, files(f).name), ...
                    'b_val', 'alpha', 'beta', 'theta_hat');
            catch
                continue;
            end

            parameter_match = ...
                S.b_val == p.b && ...
                abs(S.alpha - p.alpha) <= tol && ...
                abs(S.beta  - p.beta)  <= tol;

            if parameter_match
                chosen_file = files(f).name;
                theta_hat = S.theta_hat(:);
                break;
            end
        end

        if isempty(chosen_file)
            warning(['No estimation result found for ', ...
                'b = %d, alpha = %.2f, beta = %.2f.'], ...
                p.b, p.alpha, p.beta);
            continue;
        end


        % -----------------------------------------------------------------
        % 2.2 Construct the B-spline basis
        % -----------------------------------------------------------------

        num_basis_function = p.b + 1;
        num_breaks = num_basis_function - k_order + 2;

        breaks = linspace( ...
            Global_Min_R, ...
            Global_Max_R, ...
            num_breaks);

        knots = augknt(breaks, k_order);


        % -----------------------------------------------------------------
        % 2.3 Construct the original-grid basis used to obtain delta_t
        % -----------------------------------------------------------------

        months = Smooth_AllR.Properties.VariableNames;
        T_len = height(Realized_Return);

        Basis_Original = cell(T_len, 1);

        for t = 1:T_len
            col_name = months{t};
            R_axis_t = Smooth_AllR.(col_name);

            Basis_Original{t} = spcol( ...
                knots, ...
                k_order, ...
                R_axis_t(:));
        end


        % -----------------------------------------------------------------
        % 2.4 Recover the period-specific normalization term delta_t
        % -----------------------------------------------------------------

        R_vec_numeric  = Realized_Return.realized_ret(:);
        Rf_vec_numeric = Risk_Free_Rate(:);

        [~, ~, delta_vec] = log_likelihood_bspline( ...
            theta_hat, ...
            R_vec_numeric, ...
            Rf_vec_numeric, ...
            Basis_Original, ...
            Smooth_AllR, ...
            Smooth_AllR_RND, ...
            months, ...
            true, ...
            p.alpha, ...
            p.beta);


        % -----------------------------------------------------------------
        % 2.5 Evaluate M_t(R) on the common R grid
        % -----------------------------------------------------------------

        B_plot = spcol(knots, k_order, R_axis);

        Spline_Sum_Plot = B_plot * theta_hat;

        % Same numerical overflow protection as the original code
        Spline_Sum_Plot = max( ...
            min(Spline_Sum_Plot, 60), ...
            -60);

        % Dimensions:
        % delta_vec              : T x 1
        % Spline_Sum_Plot'       : 1 x N_R
        % M_all                  : T x N_R
        M_all = exp(delta_vec + Spline_Sum_Plot');


        % -----------------------------------------------------------------
        % 2.6 Calculate D'(F^P(R))
        % -----------------------------------------------------------------

        x = F_P_plot;

        transformed_probability = -p.beta .* log(x);

        D_plot = exp( ...
            -(transformed_probability .^ p.alpha));

        D_prime_plot = ...
            (p.alpha .* p.beta ./ x) .* ...
            (transformed_probability .^ (p.alpha - 1)) .* ...
            D_plot;

        % Protection against possible numerical problems near the endpoints
        D_prime_plot(~isfinite(D_prime_plot)) = NaN;


        % -----------------------------------------------------------------
        % 2.7 Construct small m_t(R)
        %
        % m_t(R) = M_t(R) D'(F^P(R))
        % -----------------------------------------------------------------

        m_all = M_all .* D_prime_plot';

        % Time-series average
        m_curve_avg = mean(m_all, 1, 'omitnan');


        % -----------------------------------------------------------------
        % 2.8 Plot
        % -----------------------------------------------------------------

        legend_str = sprintf( ...
            '$\\alpha=%.2f,\\;\\beta=%.2f$', ...
            p.alpha, ...
            p.beta);

        plot( ...
            R_axis, ...
            m_curve_avg, ...
            'LineWidth', 2, ...
            'Color', colors(k, :), ...
            'DisplayName', legend_str);
    end

    xlabel('$R_{t+1}$');

    ylabel( ...
        '$\mathrm{E}[m_t(R_{t+1})]$', ...
        'Interpreter', 'latex');

    grid on;
    xlim([0.8 1.2]);

    legend( ...
        'show', ...
        'Location', 'northeast', ...
        'Box', 'off', ...
        'FontSize', 11);

    out_name = sprintf( ...
        'plot_m_PricingKernel_curve_Group_%s.png', ...
        group.suffix);

    exportgraphics( ...
        fig, ...
        fullfile(Path_Output, out_name), ...
        'Resolution', 300);

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
        
        % B-spline  
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
        
        legend_str = sprintf('$\\alpha=%.2f, \\,\\beta=%.2f$', p.alpha, p.beta);
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
    exportgraphics(fig, fullfile(Path_Output, out_name), 'Resolution', 300);
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
    
    % B-Spline Derivatives
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
    % writetable(T_save, fullfile(Path_Output, csv_name));


    % ============================================================
    %  Trend Check Table (0.8, 1.0, 1.2 + 區間內 Min/Max)
    % ============================================================
    target_R_points = [0.8, 1.0, 1.2];
    trend_data = struct();
    trend_data.Measure = measures(:); 
    
    % 1. 建立 0.8 到 1.2 的遮罩，用來抓取區間內所有數據點
    range_mask = (R_axis >= 0.8) & (R_axis <= 1.2);
    
    % 2. 計算特定點 (0.8, 1.0, 1.2) 的值
    for k = 1:length(target_R_points)
        r_val = target_R_points(k);
        col_name = sprintf('Val_at_%s', strrep(num2str(r_val), '.', '_'));
        
        vals = zeros(length(measures), 1);
        for m = 1:length(measures)
            key = measures{m};
            y_data = risk_pref.(key);
            vals(m) = interp1(R_axis, y_data, r_val, 'pchip');
        end
        trend_data.(col_name) = vals;
    end
    
    % 3. 新增：找出 0.8~1.2 區間內真正的最大與最小值 (解決微笑曲線問題)
    min_vals = zeros(length(measures), 1);
    max_vals = zeros(length(measures), 1);
    for m = 1:length(measures)
        key = measures{m};
        % 只取出 R 在 0.8 到 1.2 之間的數值陣列
        y_in_range = risk_pref.(key)(range_mask);
        min_vals(m) = min(y_in_range);
        max_vals(m) = max(y_in_range);
    end
    trend_data.Min_08_12 = min_vals;
    trend_data.Max_08_12 = max_vals;
    
    T_trend = struct2table(trend_data);
    
    % 數值格式化 (保留四位小數)
    varNamesTrend = T_trend.Properties.VariableNames;
    for v = 1:numel(varNamesTrend)
        if isnumeric(T_trend.(varNamesTrend{v}))
            T_trend.(varNamesTrend{v}) = round(T_trend.(varNamesTrend{v}), 4);
        end
    end
    
    % 儲存檔案
    trend_csv_name = sprintf('TrendCheck_b_%d_alpha_%.2f_beta_%.2f.csv', p.b, p.alpha, p.beta);
    writetable(T_trend, fullfile(Path_Output, trend_csv_name));
    fprintf('Saved updated trend check: %s\n', trend_csv_name);
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
            legend_str = sprintf('$\\alpha=%.2f, \\,\\beta=%.2f$', p.alpha, p.beta);
            plot(R_axis, Y_plot, 'LineWidth', 2, 'Color', colors(k,:), 'DisplayName', legend_str);
        end
        
        xlabel('$R_{t+1}$');
        ylabel(measure_key, 'Interpreter', 'none');
        grid on;
        xlim([0.8 1.2]);
        
        % 設定 Y 軸範圍
        switch measure_key
            case 'ARA', ylim([0, 5]);   case 'RRA', ylim([0, 5]);
            case 'AP',  ylim([4, 12]);  case 'RP',  ylim([4, 12]);
            case 'AT',  ylim([4, 12]);  case 'RT',  ylim([4, 12]);
        end
        
        legend('show', 'Location', 'northeast', 'Box', 'off', 'FontSize', 11);
        
        out_name = sprintf('Group_Risk_%s_%s.png', measure_key, group.suffix);
        exportgraphics(fig, fullfile(Path_Output, out_name), 'Resolution', 300);
        close(fig);
    end
end
fprintf('Risk plots & tables done.\n');


%% Compute PIT, Plot Histogram, and Save Statistics

clc

% 設定區間數量 (bins)
num_bins = 10;

fprintf('\nGenerating PIT Histograms and Statistics...\n');
stats_list = [];
stats_cnt = 1;

use_delta = true;
T_len = length(Realized_Return.realized_ret);
months = Smooth_AllR.Properties.VariableNames;

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
    
    Basis_Precomputed = cell(T_len, 1);
    
    % Quintic B-spline 設定
    n_degree = 5;
    k_order  = n_degree + 1;    
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
    fig = figure('Position', [50 80 450 400]);
    set(fig, 'GraphicsSmoothing', 'on');
    set(fig, 'Renderer', 'painters');
    layout = tiledlayout(1, 1, 'TileSpacing', 'Compact', 'Padding', 'None');
    nexttile;
    hold on;
    
    histogram(pit_vec, num_bins, 'Normalization', 'pdf', 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'w');
    yline(1, 'r--', 'LineWidth', 2, 'DisplayName', 'Uniform(0,1)');
    hold off;
    
    xlabel('$z_t$'); ylabel('Density');
    xlim([0 1]);
    ylim([0 1.8]);
    grid on;
    set(gca, 'LooseInset', get(gca, 'TightInset'));
    
    % 統一使用 exportgraphics 並輸出至 Path_Output
    out_name = sprintf('plot_PIT_alpha_%.2f_beta_%.2f_bins_%d.png', p.alpha, p.beta, num_bins);
    exportgraphics(fig, fullfile(Path_Output, out_name), 'Resolution', 300);
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


%% Plot Risk Preference Indices (Beamer Theme: Absolute vs Relative)

% Define Color (LaTeX Beamer Theme - Metropolis)
mRed        = '#C1403D';
mLightBlue  = '#3279a8';
mDarkBlue   = '#2c3e50';
mBackground = '#FAFAFA';

% 將指標分為絕對風險與相對風險兩組
abs_measures = {'ARA', 'AP', 'AT'};
rel_measures = {'RRA', 'RP', 'RT'};

plot_sets = {
    struct('measures', {abs_measures}, 'filename', 'Group_Risk_Absolute_Beamer.png')
    struct('measures', {rel_measures}, 'filename', 'Group_Risk_Relative_Beamer.png')
};

for s = 1:2
    current_set = plot_sets{s};
    
    % Initialize figure with theme background
    fig = figure('Position', [100, 100, 1200, 400], 'Color', mBackground);
    set(fig, 'GraphicsSmoothing', 'on', 'Renderer', 'painters');
    tiledlayout(1, 3, 'TileSpacing', 'Compact', 'Padding', 'None');
    
    for m_idx = 1:3
        measure_key = current_set.measures{m_idx};
        nexttile;
        hold on;
        
        for idx = [2, 1] 
            p = param_list{idx};
            if isempty(risk_results{idx}), continue; end
            
            Y_plot = risk_results{idx}.(measure_key);
            
            % 根據參數設定，對應扭曲與非扭曲的顏色
            if idx == 1 % Distorted
                line_color = mRed;
                legend_str = sprintf('Distorted ($\\alpha=%.2f, \\,\\beta=%.2f$)', p.alpha, p.beta);
            else        % Undistorted
                line_color = mLightBlue;
                legend_str = sprintf('Undistorted ($\\alpha=%.2f, \\,\\beta=%.2f$)', p.alpha, p.beta);
            end
            
            plot(R_axis, Y_plot, 'LineWidth', 2, 'Color', line_color, 'DisplayName', legend_str);
        end
        
        % Beamer 主題標籤設定
        title(measure_key, 'FontSize', 14, 'Color', mDarkBlue, 'Interpreter', 'none');
        xlabel('$R_{t+1}$', 'Interpreter', 'latex', 'FontSize', 14);
        
        xlim([0.8 1.2]);
        switch measure_key
            case {'ARA', 'RRA'}, ylim([0, 5]);
            case {'AP', 'RP'},  ylim([4, 12]);
            case {'AT', 'RT'},  ylim([4, 12]);
        end
        grid on;
        
        % 設定 Beamer 座標軸格式
        set(gca, ...
            'box', 'on', ...
            'FontSize', 12, ...
            'XColor', mDarkBlue, ...
            'YColor', mDarkBlue, ...
            'Color', mBackground);
        
        % 只在第三張子圖加上圖例，保持畫面乾淨
        if m_idx == 3
            legend('show', 'Location', 'northeast', 'Box', 'off', 'FontSize', 11, 'TextColor', mDarkBlue, 'Interpreter', 'latex');
        end
        hold off;
    end
    
    % 匯出圖片，包含指定的背景顏色
    out_name = fullfile(Path_Output, current_set.filename);
    exportgraphics(fig, out_name, 'BackgroundColor', mBackground, 'Resolution', 300);
    close(fig);
end
fprintf('Risk plots (Beamer layout) done.\n');