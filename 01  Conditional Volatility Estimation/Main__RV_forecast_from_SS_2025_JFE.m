clear; clc;

Path_MainFolder = 'D:\Google\我的雲端硬碟\學術｜研究與論文\論文著作\MLE Pricing Kernel';
Path_Data = fullfile(Path_MainFolder, 'Code', ...
    '【2025－JFE】Conditional risk and the pricing kernel', 'data_outputs');
Path_Output = fullfile(Path_MainFolder, 'Code', '01  Output');


%% Load RV_forecast from 【2025－JFE】Conditional risk and the pricing kernel

RV_forecast_provided = load(fullfile(Path_Data, 'RV_forecast.mat'));

RV_forecast = table(RV_forecast_provided.dates_RV_forecast, RV_forecast_provided.RV_forecast_OOS, ...
    'VariableNames', {'Date', 'RV_Forecast'});


%% 

RV_forecast.Daily_Implied_MonthlyVol = ...
    (RV_forecast.RV_Forecast / sqrt(252)) * sqrt(12);
