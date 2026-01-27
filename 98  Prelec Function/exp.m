Path_MainFolder = 'D:\Google\我的雲端硬碟\學術｜研究與論文\論文著作\MLE Pricing Kernel';
Path_Output     = fullfile(Path_MainFolder, 'Code', '98  Output');

set(groot, 'defaultAxesFontName', 'Times New Roman');
set(groot, 'defaultAxesFontSize', 12);
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultLineMarkerFaceColor','auto');

x_start = -0.2;
x_end   = 1.6;
x = x_start : 0.01 : x_end;
y = exp(x);

f = figure('Position', [100, 100, 520, 500]);
layout = tiledlayout(1, 1, 'TileSpacing', 'Compact', 'Padding', 'None');
nexttile;
plot(x, y, '-', 'LineWidth', 2);
grid on;
xlabel('$x$');
ylabel('$\exp(x)$');
xlim([x_start x_end]);

FileName = sprintf('exp(x).png');
SavePath = fullfile(Path_Output, FileName);
exportgraphics(f, SavePath, 'Resolution', 300);