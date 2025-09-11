clear; clc;

Path_MainFolder = 'D:\Google\我的雲端硬碟\學術｜研究與論文\論文著作\MLE Pricing Kernel';
Path_Output = fullfile(Path_MainFolder, 'Code', '03  Output');


%% Load estimation result

mat_files = dir(fullfile(Path_Output, 'MLE_theta_b*.mat'));

for k = 1:length(mat_files)
    file_path = fullfile(Path_Output, mat_files(k).name);
    load(file_path, 'theta_hat');
    b_value = regexp(mat_files(k).name, '(?<=_b)(\d+)', 'match', 'once');
    var_name = ['theta_hat_' b_value];
    assignin('base', var_name, theta_hat);
end

clear b_value var_name k theta_hat mat_files


%% Setting

Path_Code_03 = fullfile(Path_MainFolder, 'Code', '03  MLE - Cubic B-Spline');
addpath(Path_Code_03);

R_axis = linspace(0.003, 3, 30000);

min_knot = min(R_axis);
max_knot = max(R_axis);
N = numel(R_axis);
degree = 3;                                                                % cubic B-spline


%% Plot

b_list = [4, 6, 8];

for b = b_list

    B_stack = zeros(b+1, N);
    for i = 1:(b+1)
        B_stack(i, :) = Bspline_basis_function_value( ...
            degree, b, min_knot, max_knot, i, R_axis);
    end
    
    eval(sprintf('theta_hat = theta_hat_%d;', b));
    % theta_hat = rand(1, b+1);
    
    g_vec = theta_hat * B_stack;                                           % (1×(b+1)) × ((b+1)×N) → 1×N
    
    % ===== Plot =====
    
    set(gcf, 'Position', [50, 50, 800, 400]);
    
    layout = tiledlayout(1, 2, 'TileSpacing', 'Compact', 'Padding', 'None');
    
    % ===== Left =====
    nexttile;
    hold on; grid on;
    
    plot(R_axis, B_stack, 'LineWidth', 1.5);
    plot(R_axis, g_vec, 'k', 'LineWidth', 2);
    
    hold off;
    xlabel('gross return');
    ylabel('Value');
    title('Zoomed in');
    legend([arrayfun(@(i) sprintf('$B^3_%d$', i), 0:b, 'UniformOutput', false), ...
        {'$\sum_{i=0}^{b}\theta_{i}B^3_i$'}], ...
        'Interpreter', 'latex', 'Location', 'northwest', 'Box', 'off');
    
    ylim([-0.2 2]);
    
    % ===== Right =====
    nexttile;
    hold on; grid on;
    
    plot(R_axis, B_stack, 'LineWidth', 1.5);
    plot(R_axis, g_vec, 'k', 'LineWidth', 2);
    
    hold off;
    xlabel('gross return');
    ylabel('Value');
    title('Wide range');
    % ylim([-16 3]);
    
    % ===== Save =====
    filename = sprintf('Cubic_B_Spline_b%d.png', b);
    saveas(gcf, fullfile(Path_Output, filename));

end