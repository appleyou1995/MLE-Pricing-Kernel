clear; clc; close all;

%% 1. Set paths and parameters (Same as your original settings)

Path_MainFolder = 'D:\Google\我的雲端硬碟\學術｜研究與論文\論文著作\MLE Pricing Kernel';
Path_Data = 'D:\Google\我的雲端硬碟\學術｜研究與論文\論文著作\CDI Method';

% Set TTM to check
Target_TTM = 30; 
L = 1;


%% 2. Load original data (Risk_Free_Rate & Realized_Return & RND)

fprintf('Loading original data...\n');

% Load risk-free rate
Path_Data_01 = fullfile(Path_Data, 'Code', '01  輸出資料');
FileName = 'Risk_Free_Rate.csv';
Risk_Free_Rate_All = readtable(fullfile(Path_Data_01, FileName));
switch Target_TTM
    case 30,  Risk_Free_Rate = Risk_Free_Rate_All.rf_gross_29d;
    case 60,  Risk_Free_Rate = Risk_Free_Rate_All.rf_gross_59d;
    case 90,  Risk_Free_Rate = Risk_Free_Rate_All.rf_gross_89d;
    case 180, Risk_Free_Rate = Risk_Free_Rate_All.rf_gross_179d;
    otherwise, error('No matching risk-free series');
end

% Load dates (from Realized Return file)
FileName = ['Realized_Return_TTM_', num2str(Target_TTM), '.csv'];
Realized_Return = readtable(fullfile(Path_Data_01, FileName));
T = length(Risk_Free_Rate);
dates = Realized_Return.date(1:T);
plot_dates = datetime(dates, 'ConvertFrom', 'yyyymmdd');

% Load Q-measure PDF tables
Path_Data_02 = fullfile(Path_Data, 'Code', '02  輸出資料 - no dividend'); % Note: Folder name kept as is
Smooth_AllR = [];
Smooth_AllR_RND = [];
years_to_merge = 1996:2021; 

for year = years_to_merge
    input_filename = fullfile(Path_Data_02, sprintf('TTM_%d_RND_Tables_%d.mat', Target_TTM, year));
    if exist(input_filename, 'file')
        data = load(input_filename);
        % Note: Assuming data structure matches your provision
        if isfield(data, 'Table_Smooth_AllR')
             Smooth_AllR     = [Smooth_AllR,     data.Table_Smooth_AllR];      
             Smooth_AllR_RND = [Smooth_AllR_RND, data.Table_Smooth_AllR_RND];  
        end
    end
end


%% 3. Load estimation results (gamma_hat, kappa_vec)

Path_Output = fullfile(Path_MainFolder, 'Code', '10  Output');
ResultFile = fullfile(Path_Output, sprintf('MLE_gamma_TTM_%d_L_%d.mat', Target_TTM, L));

if exist(ResultFile, 'file')
    load(ResultFile, 'gamma_hat', 'kappa_vec');
    fprintf('Estimation results loaded successfully: gamma = %.4f\n', gamma_hat(1));
else
    error('Estimation result file not found: %s. Please verify if the estimation program has been executed.', ResultFile);
end


%% 4. Back Calculation of Rf

% Rf_Implied = exp(-kappa) * E*[R^gamma]

Rf_Implied = zeros(T, 1);
gamma = gamma_hat(1);

for t = 1:T
    date_str = num2str(dates(t));
    
    % 4.1 Get current R axis and PDF
    R_axis = Smooth_AllR.(date_str);
    f_star_curve = Smooth_AllR_RND.(date_str);
    
    % 4.2 Compute polynomial part (Raw Integral)
    % When L=1, R^gamma is equivalent to exp(gamma * ln(R))
    integrand_raw = f_star_curve .* (R_axis .^ gamma); 
    integral_val_raw = trapz(R_axis, integrand_raw);
    
    % --- [新增] 加入與估計時一致的校正邏輯 ---
    % 1. 計算原始數據隱含的期望值 (Bias Check)
    EQ_R_biased = trapz(R_axis, f_star_curve .* R_axis);
    
    % 2. 計算校正因子 (Correction Ratio)
    Rf_t = Risk_Free_Rate(t);
    Correction_Ratio = Rf_t / EQ_R_biased;
    
    % 3. 應用校正 (得到 Refined Integral)
    % 這是計算 kappa 時真正用到的數值
    integral_val_refined = integral_val_raw * Correction_Ratio;
    
    % ---------------------------------------
    % 應用理論下界強制修正 (Theoretical Floor)
    if L == 1 && gamma > 1
         Theoretical_Min_Integral = Rf_t ^ gamma;
         
         % 如果數值積分小於理論值，驗證時也要用理論值
         if integral_val_refined < Theoretical_Min_Integral
             integral_val_refined = Theoretical_Min_Integral;
         end
    end
    % -----------------------

    % 4.3 Back-calculate Rf using the estimated kappa
    % kappa_t = -ln(Rf) + ln(Integral_Refined)
    % => ln(Rf) = ln(Integral_Refined) - kappa_t
    % => Rf = Integral_Refined / exp(kappa_t)
    Rf_Implied(t) = integral_val_refined / exp(kappa_vec(t));
end


%% 5. Plot Comparison

figure('Color', 'w', 'Name', 'Rf Check', 'Position', [100, 100, 800, 600]);
tiledlayout(2, 1, 'TileSpacing', 'Compact', 'Padding', 'None');

% Convert datetime to serial date numbers for robust plotting
x_num = datenum(plot_dates); 

% --- Top plot: Comparison of two lines ---
nexttile;
plot(x_num, Risk_Free_Rate, 'b-', 'LineWidth', 1.5); hold on;
plot(x_num, Rf_Implied, 'r--', 'LineWidth', 1.5);

title(sprintf('Check Consistency: True $R_f$ vs. Implied $R_f$ (TTM=%d)', Target_TTM), 'Interpreter', 'latex');
legend('True $R_f$ (Input)', 'Implied $R_f$ (Back-calculated from kappa)', 'Interpreter', 'latex', 'Box', 'off');
ylabel('Gross Risk Free Rate');
grid on;

% Set limits using numeric values
xlim([min(x_num), max(x_num)]);

% Format the x-axis ticks to show Years
datetick('x', 'yyyy', 'keeplimits'); 

% --- Bottom plot: Difference (Error) ---
nexttile;
Difference = Rf_Implied - Risk_Free_Rate;
plot(x_num, Difference, 'k');

title('Difference (Implied - True)');
ylabel('Error');
grid on;

% Set limits using numeric values
xlim([min(x_num), max(x_num)]);

% Format the x-axis ticks to show Years
datetick('x', 'yyyy', 'keeplimits');

% Display maximum error
max_err = max(abs(Difference));
fprintf('\nCheck completed.\n');
fprintf('Max absolute error: %.15f\n', max_err);

if max_err < 1e-10
    fprintf('Conclusion: The two are completely identical! Code logic and kappa definition are consistent.\n');
else
    fprintf('Conclusion: Error found. Please check data alignment or calculation details.\n');
end

saveas(gcf, fullfile(Path_Output, sprintf('Figure_Check_Rf_Consistency_TTM_%d.png', Target_TTM)));


%% 6. RND Martingale Check

% 檢查用 trapz 算出來的 E^Q[R] 是否等於 Rf
% 如果 Ratio < 1，代表積分範圍被截斷了（漏掉右尾）

Martingale_Ratio = zeros(T, 1);
for t = 1:T
    date_str = num2str(dates(t));
    R_axis = Smooth_AllR.(date_str);
    f_star = Smooth_AllR_RND.(date_str);
    
    % E^Q[R]
    EQ_R = trapz(R_axis, R_axis .* f_star);
    
    % Compare E^Q[R] and Rf
    Martingale_Ratio(t) = EQ_R / Risk_Free_Rate(t);
end

figure;
plot(plot_dates, Martingale_Ratio, 'k');
yline(1, 'r--');
title(sprintf('Martingale Test: $E^Q[R] / R_f$ (TTM=%d)', Target_TTM), 'Interpreter', 'latex');
ylabel('Ratio (Should be 1)');
ylim([0.97, 1.02]);
grid on;

saveas(gcf, fullfile(Path_Output, sprintf('Figure_Check_RND_Martingale_TTM_%d.png', Target_TTM)));