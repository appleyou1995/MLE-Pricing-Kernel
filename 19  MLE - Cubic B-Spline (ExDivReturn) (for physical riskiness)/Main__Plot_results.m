clear; clc; close all;

%% Paths and settings

Path_MainFolder = ['D:\Google\我的雲端硬碟\學術｜研究與論文\論文著作\' ...
    'MLE Pricing Kernel'];
Path_Output = fullfile(Path_MainFolder, 'Code', '19  Output - Hsieh');

TTM_List = [30, 60, 90, 180];
Plot_R = (0.8:0.001:1.2)';

set(groot, 'defaultAxesFontName', 'Times New Roman');
set(groot, 'defaultAxesFontSize', 12);
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');


%% Load compact Hsieh results

R_Grid = readtable(fullfile(Path_Output, ...
    'CubicBSpline_R_Grid.csv'));
G_Shape_All = readtable(fullfile(Path_Output, ...
    'CubicBSpline_GShape_By_TTM.csv'));
Kappa_All = readtable(fullfile(Path_Output, 'Kappa_By_TTM.csv'));

require_columns(R_Grid, {'grid_index', 'gross_return_R'}, ...
    'CubicBSpline_R_Grid.csv');
require_columns(G_Shape_All, {'TTM', 'grid_index', 'g_shape'}, ...
    'CubicBSpline_GShape_By_TTM.csv');

R_Grid = sortrows(R_Grid, 'grid_index');
Full_R = double(R_Grid.gross_return_R);
if any(~isfinite(Full_R)) || any(diff(Full_R) <= 0)
    error('The common gross-return grid must be finite and increasing.');
end
if Plot_R(1) < Full_R(1) || Plot_R(end) > Full_R(end)
    error('The requested plotting range is outside the exported R grid.');
end

M_avg_allTTM = NaN(numel(TTM_List), numel(Plot_R));


%% Reconstruct and average the monthly pricing kernels
%
% g_shape_T(R) = exp[-q_T(R)]
% M_t,T(R)      = exp(kappa_t,T + q_T(R))
%               = exp(kappa_t,T) / g_shape_T(R)

for i = 1:numel(TTM_List)
    TTM = TTM_List(i);

    Shape_Table = G_Shape_All(G_Shape_All.TTM == TTM, :);
    Shape_Table = sortrows(Shape_Table, 'grid_index');
    if height(Shape_Table) ~= height(R_Grid) || ...
            ~isequal(double(Shape_Table.grid_index), ...
            double(R_Grid.grid_index))
        error('The g-shape grid is invalid for TTM = %d.', TTM);
    end

    Full_G_Shape = double(Shape_Table.g_shape);
    if any(~isfinite(Full_G_Shape)) || any(Full_G_Shape <= 0)
        error('The g-shape values are invalid for TTM = %d.', TTM);
    end
    Plot_G_Shape = interp1(Full_R, Full_G_Shape, Plot_R, 'pchip');

    kappa_col = sprintf('kappa_TTM_%d', TTM);
    require_columns(Kappa_All, {kappa_col}, 'Kappa_By_TTM.csv');
    kappa_vec = double(Kappa_All.(kappa_col));
    kappa_vec = kappa_vec(isfinite(kappa_vec));
    if isempty(kappa_vec)
        error('No finite kappa values found for TTM = %d.', TTM);
    end

    M_avg_allTTM(i, :) = ...
        (mean(exp(kappa_vec)) ./ Plot_G_Shape)';
end


%% Figure 1: average pricing kernel

colors = lines(numel(TTM_List));
legend_labels = arrayfun(@(x) sprintf('TTM = %d', x), ...
    TTM_List, 'UniformOutput', false);

fig = figure('Color', 'w', 'Position', [200 200 600 460]);
tiledlayout(1, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = nexttile;
ax.Toolbar.Visible = 'off';
disableDefaultInteractivity(ax);
hold on; grid on; box on;

for i = 1:numel(TTM_List)
    plot(Plot_R, M_avg_allTTM(i, :), ...
        'LineWidth', 1.5, 'Color', colors(i, :));
end

xlim([Plot_R(1), Plot_R(end)]);
xlabel('gross return $R_T$');
ylabel('average pricing kernel $\bar{M}_t(R_T)$');
legend(legend_labels, 'Location', 'best', 'Box', 'off');
hold off;

Pricing_Kernel_File = fullfile(Path_Output, ...
    'Figure_Avg_Pricing_Kernel.png');
exportgraphics(fig, Pricing_Kernel_File, 'Resolution', 300);
close(fig);


%% Figure 2: kappa time series

fig = figure('Color', 'w', 'Position', [200 200 900 430]);
tiledlayout(1, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = nexttile;
ax.Toolbar.Visible = 'off';
disableDefaultInteractivity(ax);
hold on; grid on; box on;

for i = 1:numel(TTM_List)
    TTM = TTM_List(i);
    date_col = sprintf('date_TTM_%d', TTM);
    kappa_col = sprintf('kappa_TTM_%d', TTM);

    require_columns(Kappa_All, {date_col, kappa_col}, ...
        'Kappa_By_TTM.csv');

    date_values = double(Kappa_All.(date_col));
    kappa_values = double(Kappa_All.(kappa_col));
    valid = isfinite(date_values) & isfinite(kappa_values);
    plot_dates = datetime(string(date_values(valid)), ...
        'InputFormat', 'yyyyMMdd');

    plot(plot_dates, kappa_values(valid), ...
        'LineWidth', 1.3, 'Color', colors(i, :));
end

tick_years = [1996, 2000, 2005, 2010, 2015, 2020, 2024];
tick_dates = datetime(tick_years, 1, 1);
tick_dates(end) = datetime(2024, 12, 31);
xlim([datetime(1995, 10, 1), datetime(2025, 3, 31)]);
xticks(tick_dates);
xticklabels(string(tick_years));
ylabel('$\kappa_{t,T}$', 'Rotation', 0);
ylim([0.78, 1.00]);
legend(legend_labels, 'Location', 'southwest', 'Box', 'off');
hold off;

Kappa_File = fullfile(Path_Output, 'Figure_kappa.png');
exportgraphics(fig, Kappa_File, 'Resolution', 300);
close(fig);

fprintf('Saved: %s\n', Pricing_Kernel_File);
fprintf('Saved: %s\n', Kappa_File);


%% Local helper

function require_columns(Input_Table, Required_Columns, File_Label)
    missing_columns = setdiff(Required_Columns, ...
        Input_Table.Properties.VariableNames);
    if ~isempty(missing_columns)
        error('%s is missing columns: %s', File_Label, ...
            strjoin(missing_columns, ', '));
    end
end
