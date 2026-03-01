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

alpha_min  = 0.95;
alpha_max  = 1.05;
alpha_grid = alpha_min:diff:alpha_max;

beta_min  = 0.95;
beta_max  = 1.05;
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

Path_Output = fullfile(Path_MainFolder, 'Code', '13  Output');


%% Estimation: MLE over (alpha, beta, lambda) - Sextic P-Spline with CV

clc

% Add paths
Path_Code_13 = fullfile(Path_MainFolder, 'Code', ...
    '13  MLE - Sextic B-Spline with Penalty');
addpath(Path_Code_13);

% split sample
Realized_Return_front = Realized_Return(1:T1, :);
Risk_Free_Rate_front  = Risk_Free_Rate(1:T1);

% --- [Setting] P-Spline 與 Cross-Validation 參數 ---
use_delta = true;

% 固定給予充足的基底函數數量 (Sextic n=6 搭配 b=20)
b_val = 12;

% 建立 lambda 的對數網格 (例如 0.01 到 10000)
lambda_grid = logspace(-2, 4, 7);

% K摺交叉驗證 (取代文獻的 leave-2-out 以適應長時間序列)
K_folds = 3;

% 1. 建立任務列表 (Flatten loops)
tasks = [];
cnt = 1;
for a = 1:length(alpha_grid)
    for b_idx = 1:length(beta_grid)
        for l_idx = 1:length(lambda_grid)
            tasks(cnt).alpha  = alpha_grid(a);                             %#ok<SAGROW>
            tasks(cnt).beta   = beta_grid(b_idx);                          %#ok<SAGROW>
            tasks(cnt).lambda = lambda_grid(l_idx);                        %#ok<SAGROW>
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

fprintf('開始執行 %d 個任務，使用 CPU 平行運算 (含 %d-Fold CV)...\n', Total_Tasks, K_folds);

% 確保開啟 Parallel Pool
if isempty(gcp('nocreate')), parpool; end

% 3. 外層平行迴圈
parfor i = 1:Total_Tasks
    alpha  = tasks(i).alpha;
    beta   = tasks(i).beta;
    lambda = tasks(i).lambda;
    
    % 檔名加入 lambda 標籤 (使用科學記號避免檔名過長)
    outname = sprintf('MLE_PSpline_lam_%.1e_alpha_%.2f_beta_%.2f.mat', lambda, alpha, beta);
    OutputFile = fullfile(Path_Output, outname);    
    fprintf('Task Running: lam=%.1e, a=%.2f, b=%.2f\n', lambda, alpha, beta);
    
    D = Data_Const.Value;
    
    % --- Step A: 全樣本最終配適 (Full Sample Estimation) ---
    [theta_hat, log_lik, BIC] = MLE_BSpline_estimation( ...
        D.Smooth_AllR, D.Smooth_AllR_RND, ...
        D.Realized_Return, D.Risk_Free_Rate, ...
        b_val, use_delta, alpha, beta, D.Global_Min_R, D.Global_Max_R, lambda);
        
    % --- Step B: K-Fold 交叉驗證 (Cross-Validation) ---
    % 算出這組參數在樣本外的 Log-Likelihood 總和
    CV_log_lik = run_KFold_CV(D, b_val, use_delta, alpha, beta, lambda, K_folds, theta_hat);
    
    % --- Step C: 儲存結果 ---
    Result = struct();
    Result.theta_hat  = theta_hat;
    Result.log_lik    = log_lik;     % 全樣本的配適度
    Result.CV_log_lik = CV_log_lik;  % ★ 樣本外的配適度 (用於挑選最佳 lambda)
    Result.BIC        = BIC;
    Result.b_val      = b_val;
    Result.alpha      = alpha;
    Result.beta       = beta;
    Result.lambda     = lambda;
    
    parsave_result(OutputFile, Result);
end

disp('All tasks completed successfully.');


%% --- 輔助函數 1：用於 parfor 內存檔 ---

function parsave_result(fname, data_struct)
    save(fname, '-struct', 'data_struct');
end


%% --- 輔助函數 2：執行 K-Fold 交叉驗證 ---

function CV_LL = run_KFold_CV(D, b_val, use_delta, alpha, beta, lambda, K, theta_init)
    T = height(D.Realized_Return);
    
    % 固定亂數種子，確保每個任務的切分方式一致
    rng(0);
    
    % 使用 cvpartition 自動建立 K-fold 切分
    % (Statistics and Machine Learning Toolbox)
    cv = cvpartition(T, 'KFold', K);
    
    CV_LL_sum = 0;
    
    for k = 1:K
        % 直接取得第 k 摺的訓練集與測試集邏輯索引 (Logical array)
        train_idx = training(cv, k);
        test_idx  = test(cv, k);
        
        % 提取訓練集資料
        Train_Return = D.Realized_Return(train_idx, :);
        Train_Rf = D.Risk_Free_Rate(train_idx);
        months_all = D.Smooth_AllR.Properties.VariableNames;
        train_months = months_all(train_idx);
        Train_Smooth_R = D.Smooth_AllR(:, train_months);
        Train_Smooth_RND = D.Smooth_AllR_RND(:, train_months);
        
        % 提取測試集資料
        Test_Return = D.Realized_Return(test_idx, :);
        Test_Rf = D.Risk_Free_Rate(test_idx);
        test_months = months_all(test_idx);
        Test_Smooth_R = D.Smooth_AllR(:, test_months);
        Test_Smooth_RND = D.Smooth_AllR_RND(:, test_months);
        
        % 1. 在「訓練集」上配適模型 (取得受懲罰影響的 theta_hat)
        [theta_hat, ~, ~, ~, ~, ~] = MLE_BSpline_estimation(...
            Train_Smooth_R, Train_Smooth_RND, Train_Return, Train_Rf, ...
            b_val, use_delta, alpha, beta, D.Global_Min_R, D.Global_Max_R, lambda, theta_init);
            
        % 2. 在「測試集」上評估純粹的 Log-Likelihood (無 Penalty 污染)
        Test_LL = evaluate_test_LL(theta_hat, Test_Smooth_R, Test_Smooth_RND, ...
            Test_Return, Test_Rf, b_val, use_delta, alpha, beta, D.Global_Min_R, D.Global_Max_R);
            
        CV_LL_sum = CV_LL_sum + Test_LL;
    end
    CV_LL = CV_LL_sum;
end


%% --- 輔助函數 3：評估測試集的 Log-Likelihood ---

function Test_LL = evaluate_test_LL(theta, Smooth_AllR, Smooth_AllR_RND, ...
    Realized_Return, Risk_Free_Rate, ...
    b_val, use_delta, alpha, beta, Global_Min_R, Global_Max_R)

    % 動態預先計算測試集的 B-Spline 基底矩陣
    T = height(Realized_Return);
    R_vec = Realized_Return.realized_ret;
    Rf_vec = Risk_Free_Rate;
    months = Smooth_AllR.Properties.VariableNames;
    
    n_degree = 6; % Sextic B-spline
    k_order  = n_degree + 1;
    num_breaks = (b_val + 1) - k_order + 2;
    breaks = linspace(Global_Min_R, Global_Max_R, num_breaks);
    knots  = augknt(breaks, k_order);
    
    Basis_Precomputed = cell(T, 1);
    for t = 1:T
        col_name = months{t};
        R_axis = Smooth_AllR.(col_name);
        Basis_Precomputed{t} = spcol(knots, k_order, R_axis(:));
    end
    
    % 直接呼叫您原本撰寫的底層函數，算出未經懲罰的 LL
    [Test_LL, ~, ~, ~, ~] = log_likelihood_bspline(theta, R_vec, Rf_vec, ...
        Basis_Precomputed, Smooth_AllR, Smooth_AllR_RND, months, use_delta, alpha, beta);
end