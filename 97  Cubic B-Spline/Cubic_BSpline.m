%% MATLAB Code: Using augknt + spcol for Open Uniform B-spline
clear; clc; close all;

% --- 1. Setting ---
n_degree           = 3;
num_basis_function = 8;

k_order    = n_degree + 1;
num_breaks = num_basis_function - k_order + 2;
b          = num_basis_function - 1;

if num_breaks < 2
    error('Number of basis functions is too small; at least %d are required (k=%d).', k, k);
end

x_start = 0;
x_end   = 1;
x       = linspace(x_start, x_end, 1000);
breaks  = linspace(x_start, x_end, num_breaks);

% --- 2. Knots ---
knots = augknt(breaks, k_order);
fprintf('n = %d\n', n_degree);
fprintf('b = %d\n', b);
fprintf('Number of basis functions = %d\n', num_basis_function);
fprintf('Knots:\n');
disp(knots);

% --- 3. spcol ---
B          = spcol(knots, k_order, x);
num_basis  = size(B, 2); 
legend_str = cell(1, num_basis);

for i = 1:num_basis
    legend_str{i} = sprintf('$B_{%d}^{%d}(x)$', i-1, n_degree);
end

% --- 4. Plot ---
set(groot, 'defaultAxesFontName', 'Times New Roman');
set(groot, 'defaultAxesFontSize', 12);
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultLineMarkerFaceColor','auto');

figure('Position', [300, 300, 600, 400]);
plot(x, B, 'LineWidth', 2);
grid on;

xlabel('$x$');
ylabel('Value');
xlim([x_start x_end]);
ylim([0 1.2]);

legend(legend_str, 'Location', 'NorthWest', 'Box', 'off', 'NumColumns', 2);