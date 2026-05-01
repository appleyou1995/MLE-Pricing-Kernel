clear; clc;

Path_MainFolder = 'D:\Google\我的雲端硬碟\學術｜研究與論文\論文著作\MLE Pricing Kernel';
Path_Data = 'D:\Google\我的雲端硬碟\學術｜研究與論文\論文著作\CDI Method';
Path_Output = fullfile(Path_MainFolder, 'Code', '16  Output (High vs. Low Volatility) - Plot');


%% Load the data

Target_TTM = 30;

% Load realized gross returns (R_{t+1})
Path_Data = fullfile(Path_Data, 'Code', '01  輸出資料');
FileName = ['Old_Realized_Return_TTM_', num2str(Target_TTM), '.csv'];
Return = readtable(fullfile(Path_Data, FileName));

Path_Data = fullfile(Path_MainFolder, 'Code', '00  Data');
FileName = 'VIX.csv';
VIX = readtable(fullfile(Path_Data, FileName));

clear FileName Target_TTM


%% Mapping data

Return.Properties.VariableNames{1} = 'Date';
Return.Date = datetime(num2str(Return.Date), 'InputFormat', 'yyyyMMdd');

VIX.Properties.VariableNames{1} = 'Date';
VIX.Date = datetime(num2str(VIX.Date), 'InputFormat', 'yyyyMMdd');

VIX_subset = VIX(:, 1:2); 
VIX_subset.Properties.VariableNames{2} = 'VIX_Value';
[~, unique_idx] = unique(VIX_subset.Date, 'stable');
VIX_subset = VIX_subset(unique_idx, :);

VIX_Regimes = join(Return, VIX_subset, 'Keys', 'Date');
VIX_Regimes(:, 2) = [];
VIX_Regimes.Date = yyyymmdd(VIX_Regimes.Date);


%% Generate VIX_Regimes

vix_median = median(VIX_Regimes.VIX_Value, 'omitnan');
VIX_Regimes.Vol_Label = repmat({'Low'}, height(VIX_Regimes), 1);
VIX_Regimes.Vol_Label(VIX_Regimes.VIX_Value >= vix_median) = {'High'};
VIX_Regimes.Is_High_Vol = double(VIX_Regimes.VIX_Value >= vix_median);

out_csv = fullfile(Path_Data, 'VIX_Regimes.csv');
writetable(VIX_Regimes, out_csv);
fprintf('\nSaved to:\n%s\n', out_csv);


%% Plot

% setting
set(groot, 'defaultAxesFontName', 'Times New Roman');
set(groot, 'defaultAxesFontSize', 12);
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultLineMarkerFaceColor','auto');
colors = colororder;

if isnumeric(VIX_Regimes.Date)
    VIX_Regimes.Date = datetime(num2str(VIX_Regimes.Date), 'InputFormat', 'yyyyMMdd');
end

fig = figure('Position', [100, 100, 650, 400]);
layout = tiledlayout(1, 1, 'TileSpacing', 'Compact', 'Padding', 'None');
nexttile;
hold on;

x_start = min(VIX_Regimes.Date);
x_end = max(VIX_Regimes.Date);
y_max = max(VIX_Regimes.VIX_Value) * 1.1;

% 畫 VIX 主線條
plot(VIX_Regimes.Date, VIX_Regimes.VIX_Value, '-', ...
    'Color', colors(1,:), ...
    'LineWidth', 1.5);

% 畫中位數基準線
yline(vix_median, '--', ...
      'Color', colors(7,:), ...
      'LineWidth', 1.5, ...
      'Label', sprintf('Median = %.2f', vix_median), ...
      'LabelHorizontalAlignment', 'left', ...
      'LabelVerticalAlignment', 'bottom', ...
      'FontName', 'Times New Roman', ...
      'FontSize', 12);

% title('VIX Monthly Series and Volatility Regimes', 'FontWeight', 'bold');
% xlabel('Year', 'FontSize', 12);
ylabel('VIX Index', 'FontSize', 12);
xlim([x_start x_end]);
ylim([0 y_max]);
grid on;
box on;

out_name = sprintf('plot_VIX.png');
exportgraphics(fig, fullfile(Path_Output, out_name), 'Resolution', 300);
hold off;