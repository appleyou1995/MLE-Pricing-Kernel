clear; clc;

Path_MainFolder = 'D:\Google\我的雲端硬碟\學術｜研究與論文\論文著作\MLE Pricing Kernel';
Path_Data = 'D:\Google\我的雲端硬碟\學術｜研究與論文\論文著作\CDI Method';


%% Load the data

Target_TTM = 30;

% Load realized gross returns (R_{t+1})
Path_Data_01 = fullfile(Path_Data, 'Code', '01  輸出資料');
FileName = ['Old_Realized_Return_TTM_', num2str(Target_TTM), '.csv'];
Realized_Return = readtable(fullfile(Path_Data_01, FileName));

% Load risk-free rate R_f^t
Path_Data_01_main = fullfile(Path_Data, 'Code', '01  原始資料處理');
FileName = 'Risk_Free_Rate.csv';
Risk_Free_Rate_All = readtable(fullfile(Path_Data_01_main, FileName));
Risk_Free_Rate = Risk_Free_Rate_All.rf_gross_29d;

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
        Smooth_AllK = [Smooth_AllK, data.Table_Smooth_AllK];               %#ok<AGROW>
        Smooth_AllR = [Smooth_AllR, data.Table_Smooth_AllR];               %#ok<AGROW>
        Smooth_AllR_RND = [Smooth_AllR_RND, data.Table_Smooth_AllR_RND];   %#ok<AGROW>
    else
        warning('File %s does not exist.', input_filename);
    end
end

clear Path_Data_01 Path_Data_01_main Path_Data_02 Target_TTM
clear Risk_Free_Rate_All data FileName input_filename year


%% Distortion Coefficient

diff = 0.05;

alpha_min  = 0.7;
alpha_max  = 1.3;
alpha_grid = alpha_min:diff:alpha_max;

beta_min  = 0.9;
beta_max  = 1.1;
beta_grid = beta_min:diff:beta_max;


%% Split sample

Tq = width(Smooth_AllR_RND);
Realized_Return = Realized_Return(1:Tq, :);
Risk_Free_Rate  = Risk_Free_Rate(1:Tq);

T  = height(Realized_Return);

% All sample
T1 = T;
idx_valid = 1:T;

% Find the range of R_axis
Global_Min_R = 100; 
Global_Max_R = 0;
fields = fieldnames(Smooth_AllR);
for i = 1:numel(fields)
    this_field = fields{i};    
    if isempty(regexp(this_field, '^\d+$', 'once'))
        continue; 
    end    
    r_grid = Smooth_AllR.(this_field);    
    if ~isa(r_grid, 'double') || isempty(r_grid)
        continue;
    end    
    Global_Min_R = min(Global_Min_R, min(r_grid));
    Global_Max_R = max(Global_Max_R, max(r_grid));
end
Global_Min_R = Global_Min_R * 0.9; 
Global_Max_R = Global_Max_R * 1.1;

Path_Output = fullfile(Path_MainFolder, 'Code', '11  Output');


%% Stage 1: MLE over (alpha, beta, L)

clc

% Add paths
Path_Code_11 = fullfile(Path_MainFolder, 'Code', ...
    '11  MLE -  Exponential Polynomial with Distortion and Monotonicity (full range)');
addpath(Path_Code_11);

% split sample
Realized_Return_front = Realized_Return(1:T1, :);
Risk_Free_Rate_front  = Risk_Free_Rate(1:T1);

% setting
max_L     = 3;
use_delta = true;

for a = 1:length(alpha_grid)
    for b = 1:length(beta_grid)

        alpha = alpha_grid(a);
        beta  = beta_grid(b);

        for L = 1:max_L

            outname = sprintf('MLE_gamma_L_%d_alpha_%.2f_beta_%.2f.mat', L, alpha, beta);
            OutputFile = fullfile(Path_Output, outname);

            fprintf('\n--- Estimating: L = %d, alpha = %.2f, beta = %.2f ---\n', L, alpha, beta);
            t0 = tic;
            [gamma_hat, log_lik, delta_vec, M_vec] = MLE_gamma_estimation( ...
                Smooth_AllR, Smooth_AllR_RND, ...
                Realized_Return_front, Risk_Free_Rate_front, ...
                L, use_delta, alpha, beta, Global_Min_R, Global_Max_R);

            save(OutputFile, 'gamma_hat', 'log_lik', 'L', 'alpha', 'beta');

            elapsed = toc(t0);
            fprintf('Saved: %s (logLik=%.4g, %.2fs)\n', outname, log_lik, elapsed);
        end
    end
end


%% Plot setting

set(groot, 'defaultAxesFontName', 'Times New Roman');
set(groot, 'defaultAxesFontSize', 12);
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultLineMarkerFaceColor','auto');


%% Plot M curve - construct select_rows [2 choose 1]

% -------- User-specified (alpha, beta) by L --------
alpha_val = 0.95;
beta_val  = 0.90;

targets = [
    struct('L', 1, 'alpha', alpha_val, 'beta', beta_val);
    struct('L', 2, 'alpha', alpha_val, 'beta', beta_val);
    struct('L', 3, 'alpha', alpha_val, 'beta', beta_val)
];

select_rows = cell(3,1);
tol = 1e-9;
files = dir(fullfile(Path_Output, 'MLE_gamma_L_*.mat'));

for k = 1:numel(targets)
    L_i   = targets(k).L;
    a_tar = targets(k).alpha;
    b_tar = targets(k).beta;

    chosen = [];
    chosen_file = '';

    for f = 1:numel(files)
        Sfile = fullfile(Path_Output, files(f).name);
        S     = load(Sfile, 'L','alpha','beta','gamma_hat');

        if isfield(S,'L') && isfield(S,'alpha') && isfield(S,'beta')
            if S.L==L_i && abs(S.alpha-a_tar)<=tol && abs(S.beta-b_tar)<=tol
                chosen = S;
                chosen_file = files(f).name;
                break
            end
        end
    end
    if isempty(chosen)
        error('Exact file not found.');
    end
    select_rows{L_i} = table(L_i, a_tar, b_tar, string(chosen_file), ...
        'VariableNames', {'L','alpha','beta','file'});
end


%% Plot M curve - construct M and plot

clc

% Add paths
Path_Code_11 = fullfile(Path_MainFolder, 'Code', ...
    '11  MLE -  Exponential Polynomial with Distortion and Monotonicity (full range)');
addpath(Path_Code_11);

% plot
figure('Position',[50 80 450 400]);
layout = tiledlayout(1, 1, 'TileSpacing', 'Compact', 'Padding', 'None');
nexttile;
hold on;

sp = 3;  % All sample
colors = get(groot, 'defaultAxesColorOrder');

for L_i = 1:1
    row = select_rows{L_i};

    Sfile = fullfile(Path_Output, row.file);
    data  = load(Sfile, 'gamma_hat', 'L', 'alpha', 'beta');
    fprintf('\nL = %d, alpha = %.2f, beta = %.2f, gamma_hat = %s\n', ...
        row.L, row.alpha, row.beta, mat2str(data.gamma_hat, 4));

    [R_axis, M_bar] = compute_M_curve( ...
        data.gamma_hat, row.L, row.alpha, row.beta, ...
        samples{sp}.R, samples{sp}.Rf, samples{sp}.dates, ...
        Smooth_AllR, Smooth_AllR_RND);

    plot(R_axis, M_bar, 'LineWidth', 1.8, ...
        'Color', colors(L_i,:), ...
        'DisplayName', sprintf('$L=%d$', row.L));
end

hold off;
xlabel('$R$'); ylabel('$E(M)$');
legend('show','Location','northeast','Box','off');
grid on;
xlim([0.8 1.2]);
set(gca,'LooseInset',get(gca,'TightInset'));

% Output
out_png = sprintf('plot_M_curve_alpha_%.2f_beta_%.2f.png', alpha_val, beta_val);
saveas(gcf, fullfile(Path_Output, out_png));


%% Compute risk preference indices and plot (given L, alpha, beta)

clc
Path_Output_Risk = fullfile(Path_MainFolder, 'Code', '11  Output - Risk Peference');

% -------- Parameter sets (six cases) --------
param_list = {
    struct('L',1,'alpha',0.95,'beta',0.90)
    struct('L',2,'alpha',1.00,'beta',0.90)
    struct('L',3,'alpha',1.00,'beta',0.90)
    struct('L',1,'alpha',1.00,'beta',1.00)
    struct('L',2,'alpha',1.00,'beta',1.00)
    struct('L',3,'alpha',1.00,'beta',1.00)
};

% -------- Prepare sample (All sample) --------
R_all     = Realized_Return;
Rf_all    = Risk_Free_Rate(:);
dates_all = Realized_Return.date;

% -------- Storage for all results --------
all_results = cell(length(param_list),1);
measures = {'ARA','RRA','AP','RP','AT','RT'};

ymins = containers.Map();
ymaxs = containers.Map();

for k = 1:length(measures)
    ymins(measures{k}) = +inf;
    ymaxs(measures{k}) = -inf;
end

x_min = 0.8;
x_max = 1.2;


% ============================================================
%  Step 1: Compute M and RiskPref index

files = dir(fullfile(Path_Output, 'MLE_gamma_L_*.mat'));

for idx = 1:length(param_list)

    L_target     = param_list{idx}.L;
    alpha_target = param_list{idx}.alpha;
    beta_target  = param_list{idx}.beta;

    % ----- Find file -----
    chosen_file = '';
    for f = 1:numel(files)
        Sfile = fullfile(Path_Output, files(f).name);
        S     = load(Sfile, 'L','alpha','beta','gamma_hat');

        if S.L == L_target && S.alpha == alpha_target && S.beta == beta_target
            chosen_file = files(f).name;
            chosen      = S;
            break;
        end
    end

    if isempty(chosen_file)
        error('MLE file not found for L=%d alpha=%.2f beta=%.2f', ...
              L_target, alpha_target, beta_target);
    end

    fprintf('Loaded %s\n', chosen_file);


    % ----- Compute M curve -----
    [R_axis, M_bar] = compute_M_curve( ...
        chosen.gamma_hat, L_target, alpha_target, beta_target, ...
        R_all, Rf_all, dates_all, Smooth_AllR, Smooth_AllR_RND);

    % ----- Compute six risk preference indices -----
    [risk_pref, deriv] = compute_risk_pref_from_M(R_axis, M_bar);

    % ----- 只在 R ∈ [0.8, 1.2] 的區間內收集 y-limits -----
    R_axis_col = risk_pref.R_axis(:);
    in_range   = (R_axis_col >= x_min) & (R_axis_col <= x_max);

    % ----- Find y-limits -----
    for m = 1:length(measures)
        key = measures{m};
        vec = risk_pref.(key);
        vec = vec(:);
        vec_sub = vec(in_range);

        if ~isempty(vec_sub)
            ymin_curr = min(vec_sub,[],'omitnan');
            ymax_curr = max(vec_sub,[],'omitnan');

            if ~isnan(ymin_curr)
                ymins(key) = min(ymins(key), ymin_curr);
            end
            if ~isnan(ymax_curr)
                ymaxs(key) = max(ymaxs(key), ymax_curr);
            end
        end
    end

    % ----- Keep -----
    all_results{idx} = struct( ...
        'L',L_target, 'alpha',alpha_target, 'beta',beta_target, ...
        'R_axis',R_axis, 'M_bar',M_bar, ...
        'risk_pref',risk_pref, 'deriv',deriv );
end


% ============================================================
%  Step 2: Plot

colors = get(groot, 'defaultAxesColorOrder');

for idx = 1:length(param_list)

    res = all_results{idx};
    L_target     = res.L;
    alpha_target = res.alpha;
    beta_target  = res.beta;

    alpha_str = sprintf('%.2f', alpha_target);
    beta_str  = sprintf('%.2f', beta_target);

    for m = 1:length(measures)

        key    = measures{m};
        vec    = res.risk_pref.(key);
        R_axis = res.R_axis;

        mask = (R_axis >= x_min) & (R_axis <= x_max);
        R_plot = R_axis(mask);
        Y_plot = vec(mask);

        base_ymin = max(ymins(key), -2);
        base_ymax = min(ymaxs(key), 15);
        ylo = 0.95 * base_ymin;
        yhi = 1.05 * base_ymax;

        % 如果剛好 ylo == yhi（很罕見），給一點 buffer
        if ylo == yhi
            ylo = ylo - 1e-6;
            yhi = yhi + 1e-6;
        end

        fig = figure('Position',[50 80 450 400]);
        layout = tiledlayout(1, 1, 'TileSpacing', 'Compact', 'Padding', 'None');
        nexttile;
        hold on;
        plot(R_plot, Y_plot, 'LineWidth', 1.5, ...
            'Color', colors(L_target,:));
        grid on;
        hold off;

        xlabel('$R$');
        ylabel(key);

        xlim([x_min x_max]);
        ylim([ylo yhi]);

        out_png = sprintf('%s_L_%d_alpha_%s_beta_%s.png', ...
                          key, L_target, alpha_str, beta_str);

        saveas(fig, fullfile(Path_Output_Risk, out_png));
        close(fig);
    end

end


%% Debug Table

debug_L     = 3;
debug_alpha = 1.00;
debug_beta  = 0.90;

% Generate debug table
debug = make_debug_table(all_results, debug_L, debug_alpha, debug_beta, 1.10, 1.20);

% Build filename
alpha_str = sprintf('%.2f', debug_alpha);
beta_str  = sprintf('%.2f', debug_beta);

outname = sprintf('debug_L_%d_alpha_%s_beta_%s.csv', ...
                  debug_L, alpha_str, beta_str);

outfile = fullfile(Path_Output_Risk, outname);

% Round all numeric variables to 4 decimals
debug_rounded = debug;

vars = debug.Properties.VariableNames;
for i = 1:numel(vars)
    col = debug.(vars{i});
    if isnumeric(col)
        debug_rounded.(vars{i}) = round(col, 4);
    end
end

% Save CSV
writetable(debug_rounded, outfile);

fprintf('\nSaved rounded debug CSV:\n%s\n', outfile);
