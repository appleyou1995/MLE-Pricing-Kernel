clear; clc;

Path_MainFolder = 'D:\Google\我的雲端硬碟\學術｜研究與論文\論文著作\MLE Pricing Kernel';
Path_Data = 'D:\Google\我的雲端硬碟\學術｜研究與論文\論文著作\CDI Method';


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
Smooth_AllK = [];
Smooth_AllR = [];
Smooth_AllR_RND = [];

years_to_merge = 1996:2021;
for year = years_to_merge
    input_filename = fullfile(Path_Data_02, sprintf('TTM_%d_RND_Tables_%d.mat', Target_TTM, year));
    if exist(input_filename, 'file')
        data = load(input_filename);
        Smooth_AllK = [Smooth_AllK, data.Table_Smooth_AllK];               %#ok<AGROW>
        Smooth_AllR = [Smooth_AllR, data.Table_Smooth_AllR];               %#ok<AGROW>
        Smooth_AllR_RND = [Smooth_AllR_RND, data.Table_Smooth_AllR_RND];   %#ok<AGROW>
    else
        warning('File %s does not exist.', input_filename);
    end
end

clear Path_Data_01 Path_Data_01_main Path_Data_02 Target_TTM
clear Risk_Free_Rate_All data FileName input_filename year


%% Distortion Coefficient

diff = 0.05;

alpha_min  = 0.7;
alpha_max  = 1.3;
alpha_grid = alpha_min:diff:alpha_max;

beta_min  = 0.7;
beta_max  = 1.1;
beta_grid = beta_min:diff:beta_max;


%% Split sample

Tq = width(Smooth_AllR_RND);
Realized_Return = Realized_Return(1:Tq, :);
Risk_Free_Rate  = Risk_Free_Rate(1:Tq);

T  = height(Realized_Return);

% All sample
T1 = T;
idx_valid = 1:T;

% Find the range of R_axis
Global_Min_R = 100; 
Global_Max_R = 0;
fields = fieldnames(Smooth_AllR);
for i = 1:numel(fields)
    this_field = fields{i};    
    if isempty(regexp(this_field, '^\d+$', 'once'))
        continue; 
    end    
    r_grid = Smooth_AllR.(this_field);    
    if ~isa(r_grid, 'double') || isempty(r_grid)
        continue;
    end    
    Global_Min_R = min(Global_Min_R, min(r_grid));
    Global_Max_R = max(Global_Max_R, max(r_grid));
end
Global_Min_R = Global_Min_R * 0.9; 
Global_Max_R = Global_Max_R * 1.1;

Path_Output = fullfile(Path_MainFolder, 'Code', '12  Output');


%% Estimation: MLE over (alpha, beta, b) - Cubic B-Spline

clc

% Add paths
Path_Code_12 = fullfile(Path_MainFolder, 'Code', ...
    '12  MLE - Cubic B-Spline  with Distortion and Monotonicity (full range)');
addpath(Path_Code_12);

% split sample
Realized_Return_front = Realized_Return(1:T1, :);
Risk_Free_Rate_front  = Risk_Free_Rate(1:T1);

% setting
b_grid = [4, 6, 8];
use_delta = true;

% 1. 建立任務列表 (Flatten loops)
tasks = [];
cnt = 1;
for a = 1:length(alpha_grid)
    for b_idx = 1:length(beta_grid)
        for k = 1:length(b_grid)
            tasks(cnt).alpha = alpha_grid(a);
            tasks(cnt).beta  = beta_grid(b_idx);
            tasks(cnt).b     = b_grid(k);
            cnt = cnt + 1;
        end
    end
end

Total_Tasks = length(tasks);


% 2. [關鍵優化] 建立 Constant 資料物件
Data_Const = parallel.pool.Constant(struct(...
    'Smooth_AllR', Smooth_AllR, ...
    'Smooth_AllR_RND', Smooth_AllR_RND, ...
    'Realized_Return', Realized_Return_front, ...
    'Risk_Free_Rate', Risk_Free_Rate_front, ...
    'Global_Min_R', Global_Min_R, ...
    'Global_Max_R', Global_Max_R));

fprintf('開始執行 %d 個任務，使用 CPU 全速平行運算...\n', Total_Tasks);

% 確保開啟 Parallel Pool
if isempty(gcp('nocreate')), parpool; end


% 3. 外層平行迴圈
parfor i = 1:Total_Tasks
    
    % 取出任務參數
    b_val = tasks(i).b;
    alpha = tasks(i).alpha;
    beta  = tasks(i).beta;
    
    outname = sprintf('MLE_BSpline_b_%d_alpha_%.2f_beta_%.2f.mat', b_val, alpha, beta);
    OutputFile = fullfile(Path_Output, outname);    
    fprintf('Task Running: b=%d, a=%.2f, b=%.2f\n', b_val, alpha, beta);
    
    D = Data_Const.Value;
    
    t0 = tic;
    
    % 呼叫新的估計函數
    [theta_hat, log_lik, delta_vec, M_vec] = MLE_BSpline_estimation( ...
        D.Smooth_AllR, D.Smooth_AllR_RND, ...
        D.Realized_Return, D.Risk_Free_Rate, ...
        b_val, use_delta, alpha, beta, D.Global_Min_R, D.Global_Max_R);
    
    elapsed = toc(t0);
    
    % 存檔
    parsave_result(OutputFile, theta_hat, log_lik, b_val, alpha, beta, elapsed);
end

disp('所有任務完成！');

% --- 輔助函數：用於 parfor 內存檔 ---
function parsave_result(fname, theta_hat, log_lik, b_val, alpha, beta, elapsed)
    save(fname, 'theta_hat', 'log_lik', 'b_val', 'alpha', 'beta', 'elapsed');
end