clear; clc;

Path_MainFolder = 'D:\Google\我的雲端硬碟\學術｜研究與論文\論文著作\MLE Pricing Kernel';
Path_Data = 'D:\Google\我的雲端硬碟\學術｜研究與論文\論文著作\CDI Method';
Path_Output = fullfile(Path_MainFolder, 'Code', '02  Output');


%% Load the data

Target_TTM = 30;

% Load realized gross returns (R_{t+1})
Path_Data_01 = fullfile(Path_Data, 'Code', '01  輸出資料');
FileName = ['Realized_Return_TTM_', num2str(Target_TTM), '.csv'];
Realized_Return = readtable(fullfile(Path_Data_01, FileName));

% Load risk-free rate R_f^t
Path_Data_01_main = fullfile(Path_Data, 'Code', '01  原始資料處理');
FileName = 'Risk_Free_Rate.csv';
Risk_Free_Rate = readtable(fullfile(Path_Data_01_main, FileName));

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
        Smooth_AllK = [Smooth_AllK, data.Table_Smooth_AllK];
        Smooth_AllR = [Smooth_AllR, data.Table_Smooth_AllR];               % R_grid for interpolation
        Smooth_AllR_RND = [Smooth_AllR_RND, data.Table_Smooth_AllR_RND];   % f^*_t(R) on grid
    else
        warning('File %s does not exist.', input_filename);
    end
end

clear FileName input_filename year Path_Data_01 Path_Data_01_main Path_Data_02 data


%% Conditional volatility

Path_Code_01 = fullfile(Path_MainFolder, 'Code', '01  Output');
FileName = 'RV_forecast_with_MonthlySigma.csv';
RV_forecast_all = readtable(fullfile(Path_Code_01, FileName));

[~, idx_RV, ~] = intersect(RV_forecast_all.Date, Realized_Return.date);
RV_forecast = RV_forecast_all(idx_RV, [1, 3]);

clear Path_Code_01 FileName idx_RV


%% (1) Estimation - General Case: Estimate b

% Add paths
Path_Code_02 = fullfile(Path_MainFolder, 'Code', '02  Maximum Likelihood Estimation');
addpath(Path_Code_02);

% Initialize result storage
max_N = 5;
loglikelihoods_bfree = zeros(max_N, 1);
coefficients_bfree = NaN(max_N, max_N+1);

for N = 1:max_N
    fprintf('\n--- Estimating N = %d (b estimated) ---\n', N);

    estimate_b = true;
    b_fixed = NaN;

    [theta_hat, log_lik] = MLE_theta_estimation( ...
        Smooth_AllR, Smooth_AllR_RND, Realized_Return, Risk_Free_Rate, ...
        RV_forecast, N, estimate_b, b_fixed);

    loglikelihoods_bfree(N) = log_lik;
    coefficients_bfree(N, 1:N) = theta_hat(1:N)';                          % c_it
    coefficients_bfree(N, end) = theta_hat(end);                           % b
end

% Save table
output_matrix = [loglikelihoods_bfree'; coefficients_bfree'];
header = strcat('N=', string(1:max_N));
row_names_with_b = ["LogLikelihood", "c1", "c2", "c3", "c4", "c5", "b"];
result_table_bfree = array2table(round(output_matrix, 3), ...
    'VariableNames', header, ...
    'RowNames', row_names_with_b(1:size(output_matrix,1)));

writetable(result_table_bfree, fullfile(Path_Output, 'MLE_bfree.csv'), 'WriteRowNames', true);


%% (2) Estimation - Special Case: Fixed b = 0

% Add paths
Path_Code_02 = fullfile(Path_MainFolder, 'Code', '02  Maximum Likelihood Estimation');
addpath(Path_Code_02);

% Initialize result storage
max_N = 5;
loglikelihoods_b0 = zeros(max_N, 1);
coefficients_b0 = NaN(max_N);

for N = 1:max_N
    fprintf('\n--- Estimating N = %d (b = 0) ---\n', N);

    estimate_b = false;
    b_fixed = 0;

    [theta_hat, log_lik] = MLE_theta_estimation( ...
        Smooth_AllR, Smooth_AllR_RND, Realized_Return, Risk_Free_Rate, ...
        RV_forecast, N, estimate_b, b_fixed);

    loglikelihoods_b0(N) = log_lik;
    coefficients_b0(N, 1:N) = theta_hat';
end

% Save table
output_matrix = [loglikelihoods_b0'; coefficients_b0'];
header = strcat('N=', string(1:max_N));
row_names = ["LogLikelihood", "c1", "c2", "c3", "c4", "c5"];
result_table_b0 = array2table(round(output_matrix, 3), ...
    'VariableNames', header, ...
    'RowNames', row_names);

writetable(result_table_b0, fullfile(Path_Output, 'MLE_b0.csv'), 'WriteRowNames', true);


%% (3) Estimation - Special Case: Fixed b = 1

% Add paths
Path_Code_02 = fullfile(Path_MainFolder, 'Code', '02  Maximum Likelihood Estimation');
addpath(Path_Code_02);

% Initialize result storage
max_N = 5;
loglikelihoods_b1 = zeros(max_N, 1);
coefficients_b1 = NaN(max_N);

for N = 1:max_N
    fprintf('\n--- Estimating N = %d (b = 1) ---\n', N);

    estimate_b = false;
    b_fixed = 1;

    [theta_hat, log_lik] = MLE_theta_estimation( ...
        Smooth_AllR, Smooth_AllR_RND, Realized_Return, Risk_Free_Rate, ...
        RV_forecast, N, estimate_b, b_fixed);

    loglikelihoods_b1(N) = log_lik;
    coefficients_b1(N, 1:N) = theta_hat';
end

% Save table
output_matrix = [loglikelihoods_b1'; coefficients_b1'];
header = strcat('N=', string(1:max_N));
row_names = ["LogLikelihood", "c1", "c2", "c3", "c4", "c5"];
result_table_b1 = array2table(round(output_matrix, 3), ...
    'VariableNames', header, ...
    'RowNames', row_names);

writetable(result_table_b1, fullfile(Path_Output, 'MLE_b1.csv'), 'WriteRowNames', true);