clear; clc; close all;

%% Paths and settings

Path_MainFolder = ['D:\Google\我的雲端硬碟\學術｜研究與論文\論文著作\' ...
    'MLE Pricing Kernel'];
Path_Output = fullfile(Path_MainFolder, 'Code', '19  Output');

TTM_List = [30, 60, 90, 180];
R_axis = (0.8:0.001:1.2)';

set(groot, 'defaultAxesFontName', 'Times New Roman');
set(groot, 'defaultAxesFontSize', 12);
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');


%% Load compact CSV results

Theta_All = readtable(fullfile(Path_Output, 'Theta_By_TTM.csv'));
Kappa_All = readtable(fullfile(Path_Output, 'Kappa_By_TTM.csv'));

required_theta_columns = {'TTM', 'theta'};
if ~all(ismember(required_theta_columns, Theta_All.Properties.VariableNames))
    error('Theta_By_TTM.csv must contain the columns TTM and theta.');
end

M_avg_allTTM = NaN(numel(TTM_List), numel(R_axis));


%% Reconstruct and average the monthly pricing kernels

for i = 1:numel(TTM_List)
    TTM = TTM_List(i);

    theta_hat = Theta_All.theta(Theta_All.TTM == TTM);
    if isempty(theta_hat)
        error('No theta estimates found for TTM = %d.', TTM);
    end
    theta_hat = theta_hat(:);

    knot_file = fullfile(Path_Output, sprintf( ...
        'CubicBSpline_TTM%03d_Knots.csv', TTM));
    if ~isfile(knot_file)
        error('Knot file not found: %s', knot_file);
    end
    Knot_Table = readtable(knot_file);
    Knot_Table = sortrows(Knot_Table, 'knot_index');
    knots = Knot_Table.knot_value(:)';

    spline_order = numel(knots) - numel(theta_hat);
    if spline_order ~= 4
        error(['TTM = %d implies spline order %d. ' ...
            'A cubic B-spline requires order 4.'], TTM, spline_order);
    end

    B_mat = spcol(knots, spline_order, R_axis);
    spline_sum = B_mat * theta_hat;
    spline_sum = max(min(spline_sum, 60), -60);

    kappa_col = sprintf('kappa_TTM_%d', TTM);
    if ~ismember(kappa_col, Kappa_All.Properties.VariableNames)
        error('Kappa_By_TTM.csv is missing column %s.', kappa_col);
    end
    kappa_vec = double(Kappa_All.(kappa_col));
    kappa_vec = kappa_vec(isfinite(kappa_vec));
    if isempty(kappa_vec)
        error('No finite kappa values found for TTM = %d.', TTM);
    end

    log_M_all = kappa_vec + spline_sum';
    M_all = exp(log_M_all);
    M_avg_allTTM(i, :) = mean(M_all, 1);
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
    plot(R_axis, M_avg_allTTM(i, :), ...
        'LineWidth', 1.5, 'Color', colors(i, :));
end

xlim([R_axis(1), R_axis(end)]);
xlabel('gross return $R_T$');
ylabel('average pricing kernel $\bar{M}_t(R_T)$');
legend(legend_labels, 'Location', 'best', 'Box', 'off');
hold off;

exportgraphics(fig, fullfile(Path_Output, ...
    'Figure_Avg_Pricing_Kernel.png'), 'Resolution', 300);
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

    if ~all(ismember({date_col, kappa_col}, ...
            Kappa_All.Properties.VariableNames))
        error('Kappa_By_TTM.csv is missing columns for TTM = %d.', TTM);
    end

    date_values = double(Kappa_All.(date_col));
    kappa_values = double(Kappa_All.(kappa_col));
    valid = isfinite(date_values) & isfinite(kappa_values);
    plot_dates = datetime(string(date_values(valid)), ...
        'InputFormat', 'yyyyMMdd');

    plot(plot_dates, kappa_values(valid), ...
        'LineWidth', 1.3, 'Color', colors(i, :));
end

ylabel('$\kappa_t$', 'Rotation', 0);
legend(legend_labels, 'Location', 'best', 'Box', 'off');
hold off;

exportgraphics(fig, fullfile(Path_Output, ...
    'Figure_kappa.png'), 'Resolution', 300);
close(fig);

fprintf('Saved: %s\n', fullfile(Path_Output, ...
    'Figure_Avg_Pricing_Kernel.png'));
fprintf('Saved: %s\n', fullfile(Path_Output, 'Figure_kappa.png'));
