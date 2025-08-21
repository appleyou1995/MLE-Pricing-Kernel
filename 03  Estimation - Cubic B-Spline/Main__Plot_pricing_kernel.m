clear; clc;

Path_MainFolder = 'D:\Google\我的雲端硬碟\學術｜研究與論文\論文著作\MLE Pricing Kernel';
Path_Data = 'D:\Google\我的雲端硬碟\學術｜研究與論文\論文著作\CDI Method';
Path_Output = fullfile(Path_MainFolder, 'Code', '03  Output');


%% Load the data

Target_TTM = 30;

% Load realized gross returns (R_{t+1})
Path_Data_01 = fullfile(Path_Data, 'Code', '01  輸出資料');
FileName = ['Realized_Return_TTM_', num2str(Target_TTM), '.csv'];
Realized_Return = readtable(fullfile(Path_Data_01, FileName));

% Load risk-free rate R_f^t
Path_Data_01_main = fullfile(Path_Data, 'Code', '01  原始資料處理');
FileName = 'Risk_Free_Rate.csv';
Risk_Free_Rate = readtable(fullfile(Path_Data_01_main, FileName));

% Load Q-measure PDF tables: R axis and corresponding f^*_t(R)
Path_Data_02 = fullfile(Path_Data, 'Code', '02  輸出資料');
Smooth_AllK = [];
Smooth_AllR = [];
Smooth_AllR_RND = [];

years_to_merge = 1996:2021;
for year = years_to_merge
    input_filename = fullfile(Path_Data_02, sprintf('TTM_%d_RND_Tables_%d.mat', Target_TTM, year));
    if exist(input_filename, 'file')
        data = load(input_filename);
        Smooth_AllK = [Smooth_AllK, data.Table_Smooth_AllK];
        Smooth_AllR = [Smooth_AllR, data.Table_Smooth_AllR];               % R_grid for interpolation
        Smooth_AllR_RND = [Smooth_AllR_RND, data.Table_Smooth_AllR_RND];   % f^*_t(R) on grid
    else
        warning('File %s does not exist.', input_filename);
    end
end

clear FileName input_filename year Path_Data_01 Path_Data_01_main Path_Data_02 data years_to_merge


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


%% Function: compute_M_spline_on_date

function M_plot = compute_M_spline_on_date(theta, b, R_axis_full, f_star_full, R_plot, Rf_t)

    degree = 3;                                                            % cubic B-spline

    % --- Full grid basis for integral ---
    min_knot = min(R_axis_full);
    max_knot = max(R_axis_full);
    N_full   = numel(R_axis_full);

    B_full = zeros(b+1, N_full);
    for i = 1:(b+1)
        B_full(i, :) = Bspline_basis_function_value(degree, b, min_knot, max_knot, i, R_axis_full);
    end
    g_vec_full = theta(:).' * B_full;                                      % 1×N_full

    % δ_t
    integrand   = f_star_full .* exp(-g_vec_full);
    integralVal = trapz(R_axis_full, integrand);
    delta_t     = -log(Rf_t) + log(integralVal);

    % --- Plot grid basis for drawing ---
    B_plot = zeros(b+1, numel(R_plot));
    for i = 1:(b+1)
        B_plot(i, :) = Bspline_basis_function_value(degree, b, min_knot, max_knot, i, R_plot);
    end
    g_vec_plot = theta(:).' * B_plot;                                      % 1×N_plot

    % M(R) on plot grid
    M_plot = exp(delta_t + g_vec_plot);
end


%% Setting

b_list = [4, 6, 8];

% Choose a plotting time index t
t = 100;
date_str = num2str(Realized_Return.date(t));

% Pull inputs for date t
R_axis_full = Smooth_AllR.(date_str);                                      % 1×N grid of gross returns
f_star_full = Smooth_AllR_RND.(date_str);                                  % 1×N Q-measure PDF on the grid

% Risk-free R^f_t
Rf_t_annual = Risk_Free_Rate.rate(t);
Rf_t        = exp(Rf_t_annual * (Target_TTM / 252));

% Range of R_axis
idx_plot  = (R_axis_full >= 0.8) & (R_axis_full <= 1.2);
R_plot    = R_axis_full(idx_plot);
f_star_pl = f_star_full(idx_plot);


%% Build pricing kernel

% Add paths
Path_Code_03 = fullfile(Path_MainFolder, 'Code', '03  Estimation - Cubic B-Spline');
addpath(Path_Code_03);

M_stack = zeros(numel(R_plot), numel(b_list));

for j = 1:numel(b_list)

    b = b_list(j);

    theta_varname = sprintf('theta_hat_%d', b);
    theta_hat = eval(theta_varname);
    theta_hat = theta_hat(:);

    M_tmp = compute_M_spline_on_date(theta_hat, b, R_axis_full, f_star_full, R_plot, Rf_t);
    M_stack(:, j) = M_tmp(:);

end


%% Plot

figure; 
plot(R_plot, M_stack, 'LineWidth', 1.6); grid on;
legend(arrayfun(@(x) sprintf('b = %d', x), b_list, 'UniformOutput', false), ...
       'Location', 'northwest', 'Box', 'off');

xlabel('$R_{t+1}$', 'Interpreter', 'latex');
ylabel('$M(R_{t+1}; \theta)$', 'Interpreter', 'latex');

title(sprintf('Pricing Kernel (Spline) on %s', date_str), 'Interpreter', 'latex');

out_png = sprintf('Pricing_Kernels_Spline_%s.png', date_str);
saveas(gcf, fullfile(Path_Output, out_png));
disp(['Saved figure: ', fullfile(Path_Output, out_png)]);