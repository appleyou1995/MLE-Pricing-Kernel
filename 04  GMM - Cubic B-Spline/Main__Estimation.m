clear; clc;

Path_MainFolder = 'D:\Google\我的雲端硬碟\學術｜研究與論文\論文著作\MLE Pricing Kernel';
Path_Data = 'D:\Google\我的雲端硬碟\學術｜研究與論文\論文著作\CDI Method';
Path_Output = fullfile(Path_MainFolder, 'Code', '04  Output');


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


%% Define the knots for the B-spline

Aggregate_Smooth_AllR = Smooth_AllR.Variables;
ret_size = size(Smooth_AllK, 2);

% Find the minimum value for which the estimated risk-neutral densities have positive support
min_knot = min(Aggregate_Smooth_AllR);

% Find the maximum realized return within the sample
max_knot = 3;

clear Aggregate_Smooth_AllR


%% Estimation

% Add paths
Path_Code_04 = fullfile(Path_MainFolder, 'Code', '04  GMM - Cubic B-Spline');
addpath(Path_Code_04);

b_list = [4, 6, 8];

for b = b_list
    fprintf('\n=== Estimating GMM for b = %d ===\n', b);

    % Run estimation
    [theta_hat, Jval, delta_vec, num_vec, den_vec] = GMM_theta_estimation( ...
        Smooth_AllR, Smooth_AllR_RND, Realized_Return, Risk_Free_Rate, ...
        b, min_knot, max_knot);

    % Save result
    MatName = sprintf('GMM_theta_b%d.mat', b);
    save(fullfile(Path_Output, MatName), ...
        'theta_hat', 'Jval', 'delta_vec', 'num_vec', 'den_vec');

    % Console summary
    fprintf('b = %d | J(theta) = %.6f\n', b, Jval);
    disp('theta_hat ='); disp(theta_hat(:)');
end