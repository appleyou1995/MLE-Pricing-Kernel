clear; clc;

Path_MainFolder = 'D:\Google\我的雲端硬碟\學術｜研究與論文\論文著作\MLE Pricing Kernel';
Path_Data = 'D:\Google\我的雲端硬碟\學術｜研究與論文\論文著作\CDI Method';
Path_Function  = fullfile(Path_MainFolder, 'Code', '99  Function');
Path_Output = fullfile(Path_MainFolder, 'Code', '05  Output');


%% Load the data

addpath(Path_Function);
[Risk_Free_Rate, Smooth_AllR, Smooth_AllR_RND] = load_general_data;


%% Load estimation result

mat_files = dir(fullfile(Path_Output, 'MLE_theta_b*.mat'));

for k = 1:length(mat_files)
    file_path = fullfile(Path_Output, mat_files(k).name);
    load(file_path, 'M_vec');
    b_value = regexp(mat_files(k).name, '(?<=_b)(\d+)', 'match', 'once');
    var_name = ['M_vec_' b_value];
    assignin('base', var_name, M_vec);
end

clear file_path b_value var_name k M_vec mat_files


%% Print table

addpath(Path_Function);
build_MLE_summary(Path_Output);


%% Setting

b_list = [4, 6, 8];
M_all = {M_vec_4, M_vec_6, M_vec_8};
date_fields = Smooth_AllR.Properties.VariableNames;
R_base = Smooth_AllR.(date_fields{1});
T = numel(date_fields);
N = numel(R_base);


%% Plot setting

set(groot, 'defaultAxesFontName', 'Times New Roman');
set(groot, 'defaultAxesFontSize', 12);
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');


%% Plot Average Pricing Kernel (Full Range)

addpath(Path_Function);

plot_Average_Pricing_Kernel_Full_Range_3_subplot( ...
    Smooth_AllR, M_all, b_list, R_base, ...
    'XLim', [0 3], ...
    'YLim', [], ...
    'Path_Output', Path_Output);


%% Plot Average Pricing Kernel

addpath(Path_Function);

plot_Average_Pricing_Kernel(Smooth_AllR, M_all, b_list, R_base, ...
    'XLim', [0.8 1.2], ...
    'YLim', [0.8 2.5], ...
    'Path_Output', Path_Output);


%% Plot Average Pricing Kernel with Interval

addpath(Path_Function);

plot_Average_Pricing_Kernel_Extreme_Interval_3_subplot( ...
    Smooth_AllR, M_all, b_list, R_base, ...
    'XLim', [0.8 1.2], ...
    'YLim', [], ...
    'LowerPct', 1, ...
    'UpperPct', 99, ...
    'Path_Output', Path_Output);


%% Plot Average Pricing Kernel & Risk-free Rate

addpath(Path_Function);

[InvE_M_all, bad_count, date_dt] = plot_InvE_M_vs_Rf( ...
    Smooth_AllR, Smooth_AllR_RND, Risk_Free_Rate, M_all, b_list, ...
    'YLim', [], ...
    'Path_Output', Path_Output);