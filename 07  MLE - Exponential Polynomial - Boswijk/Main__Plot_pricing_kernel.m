clear; clc;

Path_MainFolder = 'D:\Google\我的雲端硬碟\學術｜研究與論文\論文著作\MLE Pricing Kernel';
Path_Data = 'D:\Google\我的雲端硬碟\學術｜研究與論文\論文著作\CDI Method';
Path_Function = fullfile(Path_MainFolder, 'Code', '99  Function');


%% Load the data

addpath(Path_Function);
[Risk_Free_Rate, Smooth_AllR, Smooth_AllR_RND] = load_general_data;


%% [ 2 choose 1 ] Load estimation result - with delta

Path_Output = fullfile(Path_MainFolder, 'Code', '07  Output - with delta');

mat_files = dir(fullfile(Path_Output, 'MLE_gamma_max_L_*.mat'));

for k = 1:length(mat_files)
    file_path = fullfile(Path_Output, mat_files(k).name);
    load(file_path, 'M_vec', 'delta_vec');
    L_value = regexp(mat_files(k).name, '(?<=_L_)(\d+)', 'match', 'once');
    var_name = ['M_vec_' L_value];
    assignin('base', var_name, M_vec);
    assignin('base', ['delta_vec_' L_value], delta_vec);
end

M_all = {M_vec_1, M_vec_2, M_vec_3};
delta_all = {delta_vec_1, delta_vec_2, delta_vec_3};

clear file_path L_value var_name k M_vec delta_vec mat_files


%% [ 2 choose 1 ] Load estimation result - without delta

Path_Output = fullfile(Path_MainFolder, 'Code', '07  Output - without delta');

mat_files = dir(fullfile(Path_Output, 'MLE_gamma_max_L_*.mat'));

for k = 1:length(mat_files)
    file_path = fullfile(Path_Output, mat_files(k).name);
    load(file_path, 'M_vec');
    L_value = regexp(mat_files(k).name, '(?<=_L_)(\d+)', 'match', 'once');
    var_name = ['M_vec_' L_value];
    assignin('base', var_name, M_vec);
end

M_all = {M_vec_1, M_vec_2, M_vec_3};

clear file_path L_value var_name k M_vec mat_files


%% Print table

addpath(Path_Function);
Tcmd = build_MLE_summary_Boswijk(Path_Output);


%% Setting

max_L_list = [1, 2, 3];
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
    Smooth_AllR, M_all, max_L_list, R_base, ...
    'XLim', [0 3], ...
    'YLim', [0 20], ...
    'Path_Output', Path_Output);


%% Plot Average Pricing Kernel

addpath(Path_Function);

plot_Average_Pricing_Kernel(Smooth_AllR, M_all, max_L_list, R_base, ...
    'XLim', [0.8 1.2], ...
    'YLim', [0.6 2.1], ...
    'Path_Output', Path_Output);


%% Plot Average Pricing Kernel with Interval

addpath(Path_Function);

plot_Average_Pricing_Kernel_Extreme_Interval_3_subplot( ...
    Smooth_AllR, M_all, max_L_list, R_base, ...
    'XLim', [0.8 1.2], ...
    'YLim', [0.99 1.01], ...
    'LowerPct', 1, ...
    'UpperPct', 99, ...
    'Path_Output', Path_Output);


%% Plot delta_t Time Series (MLE)

addpath(Path_Function);

[fig, ~] = plot_Delta_Time_Series_MLE( ...
    date_fields, delta_all, max_L_list, ...
    'Path_Output', Path_Output);