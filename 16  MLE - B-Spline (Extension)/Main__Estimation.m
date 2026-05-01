clear; clc;

Path_Main = 'D:\Google\我的雲端硬碟\學術｜研究與論文\論文著作\MLE Pricing Kernel\Code';

Path_Data    = 'D:\Google\我的雲端硬碟\學術｜研究與論文\論文著作\CDI Method';
Path_Data_01 = fullfile(Path_Data, 'Code', '01  輸出資料');
Path_Data_02 = fullfile(Path_Data, 'Code', '02  輸出資料');


%% Master Switch

Target_TTM = 30;
% [30, 60, 90, 180]

Extension_Mode = 'Full';
% ['Full', 'Non_Recession', 'High_Vol', 'Low_Vol']

fprintf('Current Mode: %s (TTM = %d)\n', Extension_Mode, Target_TTM);


%% Load Data

% Load risk-free rate (R_f^t)
FileName = 'Risk_Free_Rate.csv';
Risk_Free_Rate_All = readtable(fullfile(Path_Data_01, FileName));

switch Target_TTM
    case 30,  Risk_Free_Rate = Risk_Free_Rate_All.rf_gross_29d;
    case 60,  Risk_Free_Rate = Risk_Free_Rate_All.rf_gross_59d;
    case 90,  Risk_Free_Rate = Risk_Free_Rate_All.rf_gross_89d;
    case 180, Risk_Free_Rate = Risk_Free_Rate_All.rf_gross_179d;
    otherwise, error('No matching risk-free series for Target_TTM = %d', Target_TTM);
end

% Load realized gross returns (R_{t+1})
FileName = ['div_Realized_Return_TTM_', num2str(Target_TTM), '.csv'];
Realized_Return = readtable(fullfile(Path_Data_01, FileName));

% Load Q-measure PDF tables: R axis and corresponding f^*_t(R)
Smooth_AllK = [];
Smooth_AllR = [];
Smooth_AllR_RND = [];

for year = 1996:2021
    input_filename = fullfile(Path_Data_02, sprintf('TTM_%d_RND_Tables_%d.mat', Target_TTM, year));
    if exist(input_filename, 'file')
        data = load(input_filename);
        Smooth_AllK     = [Smooth_AllK, data.Table_Smooth_AllK];           %#ok<AGROW>
        Smooth_AllR     = [Smooth_AllR, data.Table_Smooth_AllR];           %#ok<AGROW>
        Smooth_AllR_RND = [Smooth_AllR_RND, data.Table_Smooth_AllR_RND];   %#ok<AGROW>
    else
        warning('File %s does not exist.', input_filename);
    end
end

clear data Risk_Free_Rate_All Path_Data_01 Path_Data_02 input_filename year FileName


%% Subsample Filter

months_all = Smooth_AllR.Properties.VariableNames;
T_total = length(months_all);

Realized_Return = Realized_Return(1:T_total, :);
Risk_Free_Rate  = Risk_Free_Rate(1:T_total);

idx_keep = true(T_total, 1);

switch Extension_Mode
    case 'Full'
        Path_Output = fullfile(Path_Main, '16  Output (TTM Comparison)');
        
    case 'Non_Recession'
        % NBER 官方定義的樣本內所有衰退期
        recession_yyyymm = [200103:200111, 200712, 200801:200812, 200901:200906, 202002:202004];
        
        % 從 Realized_Return.Date 提取 YYYYMM
        % current_yyyymm = year(Realized_Return.date)*100 + month(Realized_Return.date);
        current_yyyymm = floor(Realized_Return.date/100);
        idx_keep = ~ismember(current_yyyymm, recession_yyyymm);
        Path_Output = fullfile(Path_Main, '16  Output (Non-Recession Periods)');
        
    case 'High_Vol'
        VIX_Data = readtable(fullfile(Path_Main, '00  Data', 'VIX_Regimes.csv'));
        idx_keep = VIX_Data.Is_High_Vol == 1;
        Path_Output = fullfile(Path_Main, '16  Output (High vs. Low Volatility)');
        
    case 'Low_Vol'
        VIX_Data = readtable(fullfile(Path_Main, '00  Data', 'VIX_Regimes.csv'));
        idx_keep = VIX_Data.Is_High_Vol == 0;
        Path_Output = fullfile(Path_Main, '16  Output (High vs. Low Volatility)');
end

% 執行篩選
Realized_Return_Sub = Realized_Return(idx_keep, :);
Risk_Free_Rate_Sub  = Risk_Free_Rate(idx_keep);
valid_months_Sub    = months_all(idx_keep);

fprintf('樣本篩選完成，共包含 %d 個月份。\n', length(valid_months_Sub));


%% Determine Global Bounds (for Consistent B-Spline Knots)

Global_Min_R = 100;
Global_Max_R = 0;

for i = 1:numel(months_all)
    col_data = Smooth_AllR.(months_all{i}); 
    if isnumeric(col_data) && ~isempty(col_data)
        Global_Min_R = min(Global_Min_R, min(col_data));
        Global_Max_R = max(Global_Max_R, max(col_data));
    end
end

Global_Min_R = Global_Min_R * 0.9; 
Global_Max_R = Global_Max_R * 1.1;


%% Distortion Coefficient

diff = 0.05;

alpha_min  = 0.90;
alpha_max  = 1.05;
alpha_grid = alpha_min:diff:alpha_max;

beta_min  = 0.85;
beta_max  = 1.05;
beta_grid = beta_min:diff:beta_max;


%% Estimation

% Add path
Path_Code_16 = fullfile(Path_Main, '16  MLE - B-Spline (Extension)');
addpath(Path_Code_16);

% Setting
b_grid    = 6;
use_delta = true;

% Flatten loops
Total_Combinations = length(alpha_grid) * length(beta_grid) * length(b_grid);
tasks = repmat(struct('alpha', 0, 'beta', 0, 'b', 0), Total_Combinations, 1);

cnt = 1;
for a_idx = 1:length(alpha_grid)
    for b_idx = 1:length(beta_grid)
        for k = 1:length(b_grid)
            tasks(cnt).alpha = alpha_grid(a_idx);
            tasks(cnt).beta  = beta_grid(b_idx);
            tasks(cnt).b     = b_grid(k);
            cnt = cnt + 1;
        end
    end
end
Total_Tasks = length(tasks);

% 建立 Constant 資料物件
Data_Const = parallel.pool.Constant(struct(...
    'Smooth_AllR', Smooth_AllR, ...
    'Smooth_AllR_RND', Smooth_AllR_RND, ...
    'Realized_Return', Realized_Return_Sub, ...
    'Risk_Free_Rate', Risk_Free_Rate_Sub, ...
    'valid_months', {valid_months_Sub}, ...
    'Global_Min_R', Global_Min_R, ...
    'Global_Max_R', Global_Max_R));

fprintf('開始執行 %d 個任務...\n', Total_Tasks);

% 確保開啟 Parallel Pool
if isempty(gcp('nocreate')), parpool; end

% 外層平行迴圈
parfor i = 1:Total_Tasks
    b_val = tasks(i).b;
    alpha = tasks(i).alpha;
    beta  = tasks(i).beta;
    
    outname = sprintf('MLE_BSpline_TTM_%d_%s_alpha_%.2f_beta_%.2f.mat', ...
        Target_TTM, Extension_Mode, alpha, beta);
    OutputFile = fullfile(Path_Output, outname);    
    fprintf('Task Running: TTM=%d, %s, a=%.2f, b=%.2f\n', ...
        Target_TTM, Extension_Mode, alpha, beta);
    
    D = Data_Const.Value;
    
    [theta_hat, log_lik, BIC] = MLE_BSpline_estimation( ...
        D.Smooth_AllR, D.Smooth_AllR_RND, ...
        D.Realized_Return, D.Risk_Free_Rate, ...
        b_val, use_delta, alpha, beta, ...
        D.Global_Min_R, D.Global_Max_R, D.valid_months);
    
    Result = struct();
    Result.theta_hat = theta_hat;
    Result.log_lik   = log_lik;
    Result.b_val     = b_val;
    Result.alpha     = alpha;
    Result.beta      = beta;
    % Result.BIC       = BIC;
    
    parsave_result(OutputFile, Result);
end
disp('All tasks completed.');

% --- 輔助函數：用於 parfor 內存檔 ---
function parsave_result(fname, data_struct)
    save(fname, '-struct', 'data_struct');
end