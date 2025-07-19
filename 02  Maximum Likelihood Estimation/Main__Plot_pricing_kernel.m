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


%% Load estimation result

file_path = fullfile(Path_Output, 'MLE_Estimation_Results.csv');
estimation_table = readtable(file_path, 'ReadRowNames', true);

max_N = width(estimation_table);


%% Setting

R_plot = linspace(0.8, 1.2, 10000)';
logR = log(R_plot);

t = 291;

f_star_R = Smooth_AllR_RND(:, t);
f_star_R = f_star_R{1,1};

R_grid = Smooth_AllR(:, t);
R_grid = R_grid{1,1};


%% Compute delta

Rf_t = Risk_Free_Rate{t,2};  % 假設取第t期

% 準備儲存每一組 M 函數
M_all = zeros(length(R_plot), max_N);

for col = 1:max_N
    % 取得對應的係數向量
    c_names = estimation_table.Properties.RowNames(2:end);  % 去掉 LogLikelihood
    c_values = estimation_table{c_names, col};
    c_values = c_values(~isnan(c_values));  % 移除 NaN
    N = length(c_values);
    
    % 計算 inner exponential term over grid（在 R_grid 上積分）
    logR_grid = log(R_grid);
    exponent_term = zeros(size(R_grid));
    for i = 1:N
        exponent_term = exponent_term + c_values(i) .* (logR_grid.^i);
    end
    discount_factor = exp(-exponent_term);

    % 積分計算 δₜ（使用 trapezoidal rule）
    integrand = f_star_R .* discount_factor;
    integral_value = trapz(R_grid, integrand);
    delta_t = -log(Rf_t) + log(integral_value);

    % 對應每一個 R，計算 M(R)
    exponent_plot = zeros(size(R_plot));
    for i = 1:N
        exponent_plot = exponent_plot + c_values(i) .* (logR.^i);
    end
    M_all(:, col) = exp(delta_t + exponent_plot);
end


%% Plot

figure;
plot(R_plot, M_all, 'LineWidth', 1.5)
legend("N=1", "N=2", "N=3", "N=4", "N=5", 'Location', 'best')
xlabel('$R_{t+1}$', 'Interpreter', 'latex')
ylabel('$M(R_{t+1}; \theta)$', 'Interpreter', 'latex')
title('Estimated Pricing Kernels under Different $N$', 'Interpreter', 'latex')
grid on;
