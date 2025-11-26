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

% ----- All sample -----
T1 = T;
idx_valid = 1:T;
Path_Output = fullfile(Path_MainFolder, 'Code', '09  Output');


%% Stage 1: MLE over (alpha, beta, L)

clc

% Add paths
Path_Code_09 = fullfile(Path_MainFolder, 'Code', ...
    '09  MLE -  Exponential Polynomial with Distortion and Monotonicity');
addpath(Path_Code_09);

% split sample
Realized_Return_front = Realized_Return(1:T1, :);
Risk_Free_Rate_front  = Risk_Free_Rate(1:T1);

% setting
max_L     = 3;
use_delta = true;

fprintf('\n========== Stage 1: MLE (front sample only) ==========\n');
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
                L, use_delta, alpha, beta);

            save(OutputFile, 'gamma_hat', 'log_lik', 'L', 'alpha', 'beta');

            elapsed = toc(t0);
            fprintf('Saved: %s (logLik=%.4g, %.2fs)\n', outname, log_lik, elapsed);
        end
    end
end


%% Stage 2: Validation on back sample

clc

% Add paths
Path_Code_09 = fullfile(Path_MainFolder, 'Code', ...
    '09  MLE -  Exponential Polynomial with Distortion and Monotonicity');
addpath(Path_Code_09);

% split sample
Realized_Return_back = Realized_Return(idx_valid, :);
Risk_Free_Rate_back  = Risk_Free_Rate(idx_valid);
T2 = height(Realized_Return_back);

% Estimation result from stage 1
files = dir(fullfile(Path_Output, 'MLE_gamma_L_*_alpha_*.mat'));

fprintf('\n========== Stage 2: Validation (back sample) ==========\n');

Results = table('Size',[0 6], ...
    'VariableTypes', {'double','double','double','double','double','string'}, ...
    'VariableNames', {'L','alpha','beta','validation_loss','logLikStage1','file'});

best.validation_loss = inf;
best.idx             = NaN;

for ff = 1:numel(files)
    S = files(ff);
    load(fullfile(S.folder, S.name), 'gamma_hat','log_lik','L','alpha','beta');

    % if beta ~= 1
    %     continue
    % end

    % calculate validation score
    validation_loss = compute_validation_error( ...
        gamma_hat, L, Smooth_AllR, Smooth_AllR_RND, ...
        Realized_Return_back, Risk_Free_Rate_back, true, ...
        alpha, beta);
    Results = [Results; {L, alpha, beta, validation_loss, log_lik, string(S.name)}]; %#ok<AGROW>

    % find best
    if validation_loss < best.validation_loss
        best.validation_loss = validation_loss;
        best.idx   = height(Results);
    end

    fprintf('Checked: %-38s  L = %d  alpha = %.2f  beta = %.2f  validation_loss = %.4g\n', ...
        S.name, L, alpha, beta, validation_loss);
end

% output
if ~isempty(Results)
    Results = sortrows(Results, 'validation_loss');
    disp(Results(1:min(10,height(Results)), :));

    b = Results(1,:);
    fprintf('\nBest by validation:');
    fprintf('\nL=%d, alpha=%.2f, beta=%.2f, validation_loss=%.4g (MLE logLik=%.4g)\n', ...
        b.L, b.alpha, b.beta, b.validation_loss, b.logLikStage1);

    save(fullfile(Path_Output, 'GMM_validation_results.mat'), 'Results');
else
    warning('No Stage-1 result files found in %s', Path_Output);
end


%% Plot setting

set(groot, 'defaultAxesFontName', 'Times New Roman');
set(groot, 'defaultAxesFontSize', 12);
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultLineMarkerFaceColor','auto');


%% Plot log-likelihood (Stage 1)

clc
figure('Position',[50 80 600 400]);
hold on;

files = dir(fullfile(Path_Output, 'MLE_gamma_L_*.mat'));

L_values = [];

for k = 1:numel(files)
    data = load(fullfile(Path_Output, files(k).name));

    % Extract L, alpha, beta, and loglik fields
    L      = data.L;
    alpha  = data.alpha;
    beta   = data.beta;
    loglik = data.log_lik;

    % Collect data
    L_values = [L_values; L];                                              %#ok<AGROW>
    data_all(k,:) = [L, alpha, beta, loglik];
end

% Convert to table for easier manipulation
T = array2table(data_all, 'VariableNames', {'L','alpha','beta','loglik'});

% Filter beta = 1 only (if you want)
T = T(T.beta == 0.9, :);

% Plot by L
unique_L = unique(T.L(~isnan(T.L)));

for i = 1:length(unique_L)
    L_i = unique_L(i);
    idx = T.L == L_i;

    [alpha_sorted, order] = sort(T.alpha(idx));
    loglik_sorted = T.loglik(idx);
    loglik_sorted = loglik_sorted(order);

    plot(alpha_sorted, loglik_sorted, '-', 'LineWidth', 1.8, ...
        'DisplayName', sprintf('L = %d', L_i));
end

hold off;
xlabel('$\alpha$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Log-Likelihood', 'FontSize', 14);
legend('show', 'Location', 'southeast', 'Box', 'off', 'FontSize', 14);
grid on;

set(gca, 'LooseInset', get(gca, 'TightInset'));

alpha_min = min(T.alpha);
alpha_max = max(T.alpha);
xticks(linspace(alpha_min, alpha_max, 13));
xlim([alpha_min alpha_max]);

% Output
out_png = 'plot_logLikStage1_beta_0.9.png';
saveas(gcf, fullfile(Path_Output, out_png));


%% Plot log-likelihood (Stage 1) - 3D

figure('Position',[50 80 700 500]);
hold on;

unique_L = unique(T.L(~isnan(T.L)));
colors = lines(length(unique_L));

for i = 1:length(unique_L)
    L_i = unique_L(i);
    idx = T.L == L_i;
    scatter3(T.alpha(idx), T.beta(idx), T.loglik(idx), ...
        60, 'filled', 'MarkerFaceColor', colors(i,:), ...
        'DisplayName', sprintf('L = %d', L_i));
end

xlabel('$\alpha$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$\beta$', 'Interpreter', 'latex', 'FontSize', 14);
zlabel('Log-Likelihood', 'FontSize', 14);
title('Stage 1: Log-Likelihood Surface by L', 'FontSize', 14);

legend('show', 'Location', 'best', 'Box', 'off', 'FontSize', 12);
grid on;
view(45,30);
colormap('turbo');


%% Plot validation_loss (Stage 2)

figure('Position',[100 100 700 450]);
hold on;

L_values = unique(Results.L);

for i = 1:length(L_values)
    L_i = L_values(i);
    idx = Results.L == L_i;
    
    [alpha_sorted, order] = sort(Results.alpha(idx));
    val_loss_sorted = Results.validation_loss(idx);
    val_loss_sorted = val_loss_sorted(order);
    
    plot(alpha_sorted, val_loss_sorted, '--', 'LineWidth', 1.8, ...
        'DisplayName', sprintf('L = %d', L_i));
end

hold off;
xlabel('$\alpha$', 'FontSize', 14);
ylabel('Validation Loss', 'FontSize', 14);
legend('show', 'Location', 'best', 'Box', 'off', 'FontSize', 14);
grid on;

set(gca, 'LooseInset', get(gca, 'TightInset'));

alpha_min = min(Results.alpha);
alpha_max = max(Results.alpha);
xticks(linspace(alpha_min, alpha_max, 16));

% Output
out_png = sprintf('plot_validation_loss_beta_1.png');
saveas(gcf, fullfile(Path_Output, out_png));


%% Plot validation_loss (Stage 2) - 3D

figure('Position',[100 100 750 500]);
hold on;

L_values = unique(Results.L);
colors = lines(length(L_values));

for i = 1:length(L_values)
    L_i = L_values(i);
    idx = Results.L == L_i;

    scatter3(Results.alpha(idx), Results.beta(idx), Results.validation_loss(idx), ...
        70, 'filled', ...
        'MarkerFaceColor', colors(i,:), ...
        'DisplayName', sprintf('L = %d', L_i));
end

xlabel('$\alpha$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$\beta$', 'Interpreter', 'latex', 'FontSize', 14);
zlabel('Validation Loss', 'FontSize', 14);
title('Stage 2: Validation Loss Surface by L', 'FontSize', 14);

legend('show', 'Location', 'northeastoutside', 'Box', 'off', 'FontSize', 12);
grid on;
view(45, 25);

colormap('turbo');
set(gca, 'LooseInset', get(gca, 'TightInset'));


%% Plot M curve - Input of compute_M_curve

% sample split
Realized_Return_front = Realized_Return(1:T1, :);
Realized_Return_back  = Realized_Return(idx_valid, :);
Realized_Return_all   = Realized_Return;

Risk_Free_Rate_front = Risk_Free_Rate(1:T1);
Risk_Free_Rate_back  = Risk_Free_Rate(idx_valid);
Risk_Free_Rate_all   = Risk_Free_Rate(:);

dates       = Realized_Return.date;
dates_front = dates(1:T1);
dates_back  = dates(idx_valid);
dates_all   = dates;

samples = { ...
    struct('name','Front sample','R',Realized_Return_front, ...
                                 'Rf',Risk_Free_Rate_front, ...
                                 'dates',dates_front), ...
    struct('name','Back sample', 'R',Realized_Return_back , ...
                                 'Rf',Risk_Free_Rate_back , ...
                                 'dates',dates_back ), ...
    struct('name','All sample' , 'R',Realized_Return, ...
                                 'Rf',Risk_Free_Rate, ...
                                 'dates',dates_all ) ...
};


%% Plot M curve - construct select_rows [2 choose 1]

% -------- Find the smallest validation_loss --------
select_rows = cell(3,1);
for L_i = 1:3
    rows = Results(Results.L==L_i & Results.beta==1, :);
    [~, ix] = min(rows.validation_loss);
    select_rows{L_i} = rows(ix, :);
end


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
Path_Code_09 = fullfile(Path_MainFolder, 'Code', ...
    '09  MLE -  Exponential Polynomial with Distortion and Monotonicity');
addpath(Path_Code_09);

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
Path_Output_Risk = fullfile(Path_MainFolder, 'Code', '09  Output - Risk Peference');

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
