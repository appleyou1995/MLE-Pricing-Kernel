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


%% Estimate θ = (c_1, ..., c_N) by maximizing log-likelihood

% Add paths
Path_Code_02 = fullfile(Path_MainFolder, 'Code', '02  Maximum Likelihood Estimation');
addpath(Path_Code_02);

% Initialize result storage
max_N = 5;
loglikelihoods = zeros(max_N, 1);
coefficients = NaN(max_N);


for N = 1:max_N
    fprintf('\n=== Estimating MLE for N = %d ===\n', N);

    [theta_hat, log_lik] = MLE_theta_estimation( ...
        Smooth_AllR, Smooth_AllR_RND, Realized_Return, Risk_Free_Rate, N);

    log_likelihoods(N) = log_lik;
    coefficients(N, 1:N) = theta_hat';
end


%% Output table

coefficients_T = coefficients';
output_matrix = [log_likelihoods; coefficients_T];

header = strcat('N', string(1:max_N));
row_names = ["LogLikelihood", "c1", "c2", "c3", "c4", "c5"];
row_names = row_names(1:size(output_matrix,1));
result_table = array2table(output_matrix, ...
    'VariableNames', header, ...
    'RowNames', row_names);

output_file = fullfile(Path_Output, 'MLE_Estimation_Results.csv');
writetable(result_table, output_file, 'WriteRowNames', true);