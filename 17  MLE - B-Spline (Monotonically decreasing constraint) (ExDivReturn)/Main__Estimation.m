clear; clc;

Path_MainFolder = 'D:\Google\我的雲端硬碟\學術｜研究與論文\論文著作\MLE Pricing Kernel';
Path_Data   = fullfile(Path_MainFolder, 'Code', '00  Output');
Path_Output = fullfile(Path_MainFolder, 'Code', '17  Output - wide');


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
Smooth_AllK = [];
Smooth_AllR = [];
Smooth_AllR_RND = [];

years_to_merge = 1996:2025;
for year = years_to_merge
    input_filename = fullfile(Path_RND, sprintf('TTM_%d_RND_Tables_%d.mat', Target_TTM, year));
    if exist(input_filename, 'file')
        data = load(input_filename);
        Smooth_AllK = [Smooth_AllK, data.Table_Smooth_AllK];               %#ok<AGROW>
        Smooth_AllR = [Smooth_AllR, data.Table_Smooth_AllR];               %#ok<AGROW>
        Smooth_AllR_RND = [Smooth_AllR_RND, data.Table_Smooth_AllR_RND];   %#ok<AGROW>
    else
        warning('File %s does not exist.', input_filename);
    end
end

clear Path_RND
clear data FileName input_filename year


%% Distortion Coefficient

diff = 0.1;

alpha_min  = 0.4;
alpha_max  = 1.3;
alpha_grid = alpha_min:diff:alpha_max;

beta_min  = 0.4;
beta_max  = 1.3;
beta_grid = beta_min:diff:beta_max;


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
Smooth_AllK     = Smooth_AllK(:, idx_rnd);
Smooth_AllR     = Smooth_AllR(:, idx_rnd);
Smooth_AllR_RND = Smooth_AllR_RND(:, idx_rnd);
Realized_Return = Realized_Return(idx_ret, :);
Risk_Free_Rate  = rf_values_all(idx_rf);

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


%% Estimation: MLE over (alpha, beta, b) - B-Spline

clc

% Add paths
Path_Code_17 = fullfile(Path_MainFolder, 'Code', ...
    '17  MLE - B-Spline (Monotonically decreasing constraint) (ExDivReturn)');
addpath(Path_Code_17);

% setting
b_grid = 6;

% 1. 建立任務列表 (Flatten loops)
tasks = [];
cnt = 1;
for a_idx = 1:length(alpha_grid)
    for b_idx = 1:length(beta_grid)
        for k = 1:length(b_grid)
            tasks(cnt).alpha = alpha_grid(a_idx);                          %#ok<SAGROW>
            tasks(cnt).beta  = beta_grid(b_idx);                           %#ok<SAGROW>
            tasks(cnt).b     = b_grid(k);                                  %#ok<SAGROW>
            cnt = cnt + 1;
        end
    end
end

Total_Tasks = length(tasks);


% 2. [關鍵優化] 建立 Constant 資料物件
Data_Const = parallel.pool.Constant(struct(...
    'Smooth_AllR', Smooth_AllR, ...
    'Smooth_AllR_RND', Smooth_AllR_RND, ...
    'Realized_Return', Realized_Return, ...
    'Risk_Free_Rate', Risk_Free_Rate, ...
    'Global_Min_R', Global_Min_R, ...
    'Global_Max_R', Global_Max_R));

fprintf('開始執行 %d 個任務，使用 CPU 平行運算...\n', Total_Tasks);

% 確保開啟 Parallel Pool
if isempty(gcp('nocreate')), parpool; end


% 3. 外層平行迴圈
parfor i = 1:Total_Tasks
    b_val = tasks(i).b;
    alpha = tasks(i).alpha;
    beta  = tasks(i).beta;
    
    outname = sprintf('MLE_BSpline_b_%d_alpha_%.2f_beta_%.2f.mat', b_val, alpha, beta);
    OutputFile = fullfile(Path_Output, outname);    
    fprintf('Task Running: b=%d, a=%.2f, b=%.2f\n', b_val, alpha, beta);
    
    D = Data_Const.Value;
    
    [theta_hat, log_lik, BIC, exitflag, output] = MLE_BSpline_estimation( ...
        D.Smooth_AllR, D.Smooth_AllR_RND, ...
        D.Realized_Return, D.Risk_Free_Rate, ...
        b_val, alpha, beta, D.Global_Min_R, D.Global_Max_R);
    
    Result = struct();
    Result.theta_hat = theta_hat;
    Result.log_lik   = log_lik;
    Result.BIC       = BIC;
    Result.b_val     = b_val;
    Result.alpha     = alpha;
    Result.beta      = beta;
    Result.exitflag  = exitflag;
    Result.output    = output;
    
    parsave_result(OutputFile, Result);
end

disp('Done.');

% --- 輔助函數：用於 parfor 內存檔 ---
function parsave_result(fname, data_struct)
    save(fname, '-struct', 'data_struct');
end