% =========================================================================
% Prelec (1998) Probability Weighting Function Plotting
% =========================================================================

Path_MainFolder = 'D:\Google\我的雲端硬碟\學術｜研究與論文\論文著作\MLE Pricing Kernel';
Path_Output     = fullfile(Path_MainFolder, 'Code', '98  Output');

x = 0.0001:0.001:1;
alphas = [0.5, 1.0, 1.5];
betas  = [0.5, 1.0, 1.5];

set(groot, 'defaultAxesFontName', 'Times New Roman');
set(groot, 'defaultAxesFontSize', 12);
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultLineMarkerFaceColor','auto');
colors = get(groot, 'defaultAxesColorOrder');

for i = 1:length(alphas)
    for j = 1:length(betas)
        
        alpha_val = alphas(i);
        beta_val  = betas(j);
        
        % Prelec: D(x) = exp( - ( -beta * ln(x) )^alpha )
        y = exp( - ( -beta_val .* log(x) ).^alpha_val );
        
        f = figure('Position', [100, 100, 520, 500]);
        layout = tiledlayout(1, 1, 'TileSpacing', 'Compact', 'Padding', 'None');
        nexttile;
        hold on;
        plot(x, y, '-', 'LineWidth', 2, 'Color', colors(1,:)); hold on;        
        plot([0, 1], [0, 1], '--', 'LineWidth', 1.5, 'Color', colors(7,:));
        
        xlabel('$x$');
        ylabel('$D(x)$');
        
        xlim([0 1]);
        ylim([0 1]);
        grid on;
        axis square;
        
        FileName = sprintf('Prelec_Alpha_%.1f_Beta_%.1f.png', alpha_val, beta_val);
        SavePath = fullfile(Path_Output, FileName);
        
        exportgraphics(f, SavePath, 'Resolution', 300);        
        fprintf('Saved: %s\n', FileName);        
    end
end