clear; clc;

Path_MainFolder = 'D:\Google\我的雲端硬碟\學術｜研究與論文\論文著作\MLE Pricing Kernel';
Path_Data = 'D:\Google\我的雲端硬碟\學術｜研究與論文\論文著作\CDI Method';


%% Load the data

% Target_TTM = [30, 60, 90, 180]
Target_TTM = 30;

% Load risk-free rate R_f^t
Path_Data_01 = fullfile(Path_Data, 'Code', '01  輸出資料');
FileName = 'Risk_Free_Rate.csv';
Risk_Free_Rate_All = readtable(fullfile(Path_Data_01, FileName));

switch Target_TTM
    case 30
        Risk_Free_Rate = Risk_Free_Rate_All.rf_gross_29d;
    case 60
        Risk_Free_Rate = Risk_Free_Rate_All.rf_gross_59d;
    case 90
        Risk_Free_Rate = Risk_Free_Rate_All.rf_gross_89d;
    case 180
        Risk_Free_Rate = Risk_Free_Rate_All.rf_gross_179d;
    otherwise
        error('No matching risk-free series for Target_TTM = %d', Target_TTM);
end

% Load realized gross returns (R_{t+1})
Path_Data_01 = fullfile(Path_Data, 'Code', '01  輸出資料');
FileName = ['Realized_Return_TTM_', num2str(Target_TTM), '.csv'];
Realized_Return = readtable(fullfile(Path_Data_01, FileName));
T = length(Risk_Free_Rate);
Realized_Return = Realized_Return(1:T, :);

% Load Q-measure PDF tables: R axis and corresponding f^*_t(R)
Path_Data_02 = fullfile(Path_Data, 'Code', '02  輸出資料 - no dividend');
Smooth_AllK = [];
Smooth_AllR = [];
Smooth_AllR_RND = [];

years_to_merge = 1996:2021;
for year = years_to_merge
    input_filename = fullfile(Path_Data_02, sprintf('TTM_%d_RND_Tables_%d.mat', Target_TTM, year));
    if exist(input_filename, 'file')
        data = load(input_filename);
        Smooth_AllK     = [Smooth_AllK,     data.Table_Smooth_AllK];       %#ok<AGROW>
        Smooth_AllR     = [Smooth_AllR,     data.Table_Smooth_AllR];       %#ok<AGROW>
        Smooth_AllR_RND = [Smooth_AllR_RND, data.Table_Smooth_AllR_RND];   %#ok<AGROW>
    else
        warning('File %s does not exist.', input_filename);
    end
end

clear Path_Data_01 Path_Data_01_main Path_Data_02
clear Risk_Free_Rate_All years_to_merge data FileName input_filename year T


%% Estimation - with kappa

Path_Output = fullfile(Path_MainFolder, 'Code', '10  Output');

% Add paths
Path_Code_10 = fullfile(Path_MainFolder, 'Code', ...
    '10  MLE - Exponential Polynomial - Boswijk (for physical riskiness)');
addpath(Path_Code_10);

% Initialize result storage
max_L = 1;

for L = 1:max_L
    fprintf('\n--- Estimating MLE with kappa, TTM = %d ---\n', Target_TTM);

    [gamma_hat, log_lik, kappa_vec, M_cell] = MLE_gamma_estimation( ...
        Smooth_AllR, Smooth_AllR_RND, Realized_Return, Risk_Free_Rate, L);

    OutputFile = fullfile(Path_Output, sprintf('MLE_gamma_TTM_%d_L_%d.mat', Target_TTM, L));
    save(OutputFile, 'gamma_hat', 'log_lik', 'kappa_vec');
end