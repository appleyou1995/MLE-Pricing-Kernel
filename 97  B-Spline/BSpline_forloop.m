%% MATLAB Code: Loop for Open Uniform B-spline
clear; 
clc; 
close all;

% --- 0. Path Setting ---
Path_MainFolder = 'D:\Google\我的雲端硬碟\學術｜研究與論文\論文著作\MLE Pricing Kernel';
Path_Output = fullfile(Path_MainFolder, 'Code', '97  B-Spline');

% --- 1. Global Parameters ---
n_degree = 3;
k_order  = n_degree + 1;
x_start  = 0;
x_end    = 10;
x        = linspace(x_start, x_end, 1000);

switch n_degree
    case 3
        degree_name = 'Cubic';
    case 4
        degree_name = 'Quartic';
    case 5
        degree_name = 'Quintic';
    case 6
        degree_name = 'Sextic';
    otherwise
        degree_name = sprintf('Degree_%d', n_degree); 
end

% 設定繪圖風格
set(groot, 'defaultAxesFontName', 'Times New Roman');
set(groot, 'defaultAxesFontSize', 12);
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultLineMarkerFaceColor','auto');

% --- 2. Initialize Storage for Knots Table ---
target_basis_counts = (n_degree + 2) : (n_degree + 5);
num_cases = length(target_basis_counts);

% 預先計算最大 Knots 長度以建立矩陣 (Max Basis + Order)
max_knots_len = max(target_basis_counts) + k_order; 

% 建立結構陣列或矩陣來存資料
Data_List = struct();

% --- 3. Loop ---
for i = 1:num_cases
    num_basis_function = target_basis_counts(i);
    
    % Calculation
    num_breaks = num_basis_function - k_order + 2;
    b          = num_basis_function - 1; % Index usually goes 0 to b
    
    % Breaks and Knots
    breaks = linspace(x_start, x_end, num_breaks);
    knots  = augknt(breaks, k_order);
    
    % Store Knots Data (Padding with NaN for table alignment)
    current_knots_len = length(knots);
    knots_padded = nan(1, max_knots_len);
    knots_padded(1:current_knots_len) = knots;
    
    Data_List(i).Num_Basis = num_basis_function;
    Data_List(i).b_index   = b;
    Data_List(i).Knots     = knots_padded;
    
    % Basis Matrix
    B = spcol(knots, k_order, x);
    
    % Generate Legend Strings
    legend_str = cell(1, num_basis_function);
    for j = 1:num_basis_function
        legend_str{j} = sprintf('$B_{%d}^{%d}(x)$', j-1, n_degree);
    end
    
    % Plotting
    fig = figure('Position', [50 80 450 400], 'Visible', 'on');
    layout = tiledlayout(1, 1, 'TileSpacing', 'Compact', 'Padding', 'None');
    nexttile;
    plot(x, B, 'LineWidth', 2);
    grid on;
    xlabel('$x$');
    ylabel('Value');
    xlim([x_start x_end]);
    ylim([0 1.1]);
    
    % Legend settings
    if num_basis_function > 8
        n_col = 3;
    else
        n_col = 2;
    end
    legend(legend_str, 'Location', 'NorthWest', 'Box', 'off', 'NumColumns', n_col);
    
    % Save Plot
    filename_png = sprintf('%s_BSpline_b_%d.png', degree_name, b);
    full_path_png = fullfile(Path_Output, filename_png);
    saveas(fig, full_path_png);
    
    fprintf('Saved plot: %s\n', filename_png);
end

% --- 4. Generate and Save Table ---
T_info = table();
T_info.Num_Basis = [Data_List.Num_Basis]';
T_info.b_index   = [Data_List.b_index]';

Knots_Mat = vertcat(Data_List.Knots);
for k = 1:max_knots_len
    col_name = sprintf('Knot_%d', k);
    T_info.(col_name) = Knots_Mat(:, k);
end

csv_name = sprintf('%s_Knots_Summary.csv', degree_name);
table_filename = fullfile(Path_Output, csv_name);
writetable(T_info, table_filename);
fprintf('Saved table: %s\n', table_filename);

fprintf('All done.\n');