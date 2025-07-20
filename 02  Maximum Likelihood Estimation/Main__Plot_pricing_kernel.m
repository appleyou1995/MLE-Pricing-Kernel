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

t = 312;

R_grid = Smooth_AllR{1, t};
f_star_R = Smooth_AllR_RND{1, t};

r_annual = Risk_Free_Rate{t, 3};
Rf_t = exp(r_annual * (Target_TTM / 252));

plot_date = Risk_Free_Rate{t, 1};

idx = (R_grid >= 0.8) & (R_grid <= 1.2);
R_plot = R_grid(idx);
f_star_R = f_star_R(idx);

logR = log(R_plot);


%% Compute delta

M_all = zeros(length(R_plot), max_N);

for col = 1:max_N
    c_names = estimation_table.Properties.RowNames(2:end);
    c_values = estimation_table{c_names, col};
    c_values = c_values(~isnan(c_values));
    N = length(c_values);
    
    % delta_t
    logR_grid = log(R_grid);
    exponent_term = zeros(size(R_grid));
    for i = 1:N
        exponent_term = exponent_term + c_values(i) .* (logR_grid.^i);
    end
    discount_factor = exp(-exponent_term);
    integrand = Smooth_AllR_RND{1, t} .* discount_factor;
    integral_value = trapz(R_grid, integrand);
    delta_t = -log(Rf_t) + log(integral_value);

    % M(R)
    exponent_plot = zeros(size(R_plot));
    for i = 1:N
        exponent_plot = exponent_plot + c_values(i) .* (logR.^i);
    end
    M_all(:, col) = exp(delta_t + exponent_plot);
end


%% Plot

plot_date_dt = datetime(num2str(plot_date), 'InputFormat', 'yyyyMMdd');
date_str = num2str(plot_date);

figure;
plot(R_plot, M_all, 'LineWidth', 1.5)
legend("N=1", "N=2", "N=3", "N=4", "N=5", 'Location', 'northwest', 'Box', 'off')
xlabel('$R_{t+1}$', 'Interpreter', 'latex')
ylabel('$M(R_{t+1}; \theta)$', 'Interpreter', 'latex')
title(['Pricing Kernels on ', date_str, ' under Different N'], ...
      'Interpreter', 'latex')
grid on;

filename = sprintf('Pricing_Kernels_Special_Case_b=0_%s.png', date_str);
saveas(gcf, fullfile(Path_Output, filename));
clear filename