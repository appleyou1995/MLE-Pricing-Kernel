clear; clc;

Path_MainFolder = 'D:\Google\我的雲端硬碟\學術｜研究與論文\論文著作\MLE Pricing Kernel';
Path_Data = 'D:\Google\我的雲端硬碟\學術｜研究與論文\論文著作\CDI Method';


%% Load the data

Target_TTM = 30;

% Load realized gross returns (R_{t+1})
Path_Data_01 = fullfile(Path_Data, 'Code', '01  輸出資料');
FileName = ['Realized_Return_TTM_', num2str(Target_TTM), '.csv'];
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
        Smooth_AllK = [Smooth_AllK, data.Table_Smooth_AllK];
        Smooth_AllR = [Smooth_AllR, data.Table_Smooth_AllR];               % R_grid for interpolation
        Smooth_AllR_RND = [Smooth_AllR_RND, data.Table_Smooth_AllR_RND];   % f^*_t(R) on grid
    else
        warning('File %s does not exist.', input_filename);
    end
end

clear Path_Data_01 Path_Data_01_main Path_Data_02 Target_TTM
clear Risk_Free_Rate_All years_to_merge data FileName input_filename year


%% Distortion Coefficient

diff = 0.1;

alpha_min  = 0.5;
alpha_max  = 2.0;
alpha_grid = alpha_min:diff:alpha_max;

beta_min  = 1.0;
beta_max  = 1.0;
beta_grid = beta_min:diff:beta_max;


%% Split sample

T  = height(Realized_Return);
T1 = floor(T/2);

% validation
idx_valid = (T1+1):T;


%% First Stage Only: MLE over (alpha, beta, L)

clc
Path_Output = fullfile(Path_MainFolder, 'Code', '08  Output');

% Add paths
Path_Code_08 = fullfile(Path_MainFolder, 'Code', ...
    '08  MLE -  Exponential Polynomial with Distortion');
addpath(Path_Code_08);

% split sample
Realized_Return_front = Realized_Return(1:T1, :);
Risk_Free_Rate_front  = Risk_Free_Rate(1:T1);

% setting
max_L     = 3;
use_delta = true;
MLE_rows  = [];

fprintf('\n========== Stage 1: MLE (front sample only) ==========\n');
for a = 1:length(alpha_grid)
    for b = 1:length(beta_grid)

        alpha = alpha_grid(a);
        beta  = beta_grid(b);

        for L = 1:max_L

            outname = sprintf('MLE_gamma_L_%d_alpha_%.1f_beta_%.1f.mat', L, alpha, beta);
            OutputFile = fullfile(Path_Output, outname);

            fprintf('\n--- Estimating: L = %d, alpha = %.1f, beta = %.1f ---\n', L, alpha, beta);
            t0 = tic;
            [gamma_hat, log_lik, delta_vec, M_vec] = MLE_gamma_estimation( ...
                Smooth_AllR, Smooth_AllR_RND, ...
                Realized_Return_front, Risk_Free_Rate_front, ...
                L, use_delta, alpha, beta);

            save(OutputFile, 'gamma_hat', 'log_lik', 'L', 'alpha', 'beta');

            elapsed = toc(t0);
            fprintf('Saved: %s (logLik=%.4g, %.2fs)\n', outname, log_lik, elapsed);

            MLE_rows = [MLE_rows; {alpha, beta, L, log_lik, 'ok', OutputFile}];
        end
    end
end