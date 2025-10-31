clear; clc;

Path_MainFolder = 'D:\Google\我的雲端硬碟\學術｜研究與論文\論文著作\MLE Pricing Kernel';
Path_Data = 'D:\Google\我的雲端硬碟\學術｜研究與論文\論文著作\CDI Method';


%% Load the data

Target_TTM = 30;

% Load realized gross returns (R_{t+1})
Path_Data_01 = fullfile(Path_Data, 'Code', '01  輸出資料');
FileName = ['Realized_Return_TTM_', num2str(Target_TTM), '.csv'];
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
clear Risk_Free_Rate_All years_to_merge data FileName input_filename year


%% Distortion Coefficient

diff = 0.1;

alpha_min  = 0.5;
alpha_max  = 2.0;
alpha_grid = alpha_min:diff:alpha_max;

beta_min  = 0.5;
beta_max  = 0.9;
beta_grid = beta_min:diff:beta_max;


%% Split sample

T  = height(Realized_Return);
T1 = floor(T/2);

% validation
idx_valid = (T1+1):T;


%% Stage 1: MLE over (alpha, beta, L)

clc
Path_Output = fullfile(Path_MainFolder, 'Code', '08  Output');

% Add paths
Path_Code_08 = fullfile(Path_MainFolder, 'Code', ...
    '08  MLE -  Exponential Polynomial with Distortion');
addpath(Path_Code_08);

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

            outname = sprintf('MLE_gamma_L_%d_alpha_%.1f_beta_%.1f.mat', L, alpha, beta);
            OutputFile = fullfile(Path_Output, outname);

            fprintf('\n--- Estimating: L = %d, alpha = %.1f, beta = %.1f ---\n', L, alpha, beta);
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
Path_Output = fullfile(Path_MainFolder, 'Code', '08  Output');

% Add paths
Path_Code_08 = fullfile(Path_MainFolder, 'Code', ...
    '08  MLE -  Exponential Polynomial with Distortion');
addpath(Path_Code_08);

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

    if beta ~= 1
        continue
    end

    % calculate validation score
    validation_loss = Compute_validation_error( ...
        gamma_hat, L, Smooth_AllR, Smooth_AllR_RND, ...
        Realized_Return_back, Risk_Free_Rate_back, true, ...
        alpha, beta);
    Results = [Results; {L, alpha, beta, validation_loss, log_lik, string(S.name)}]; %#ok<AGROW>

    % find best
    if validation_loss < best.validation_loss
        best.validation_loss = validation_loss;
        best.idx   = height(Results);
    end

    fprintf('Checked: %-38s  L = %d  alpha = %.1f  beta = %.1f  validation_loss = %.4g\n', ...
        S.name, L, alpha, beta, validation_loss);
end

% output
if ~isempty(Results)
    Results = sortrows(Results, 'validation_loss');
    disp(Results(1:min(10,height(Results)), :));

    b = Results(1,:);
    fprintf('\nBest by validation:');
    fprintf('\nL=%d, alpha=%.1f, beta=%.1f, validation_loss=%.4g (stage1 logLik=%.4g)\n', ...
        b.L, b.alpha, b.beta, b.validation_loss, b.logLikStage1);

    save(fullfile(Path_Output, 'Stage2_validation_results.mat'), 'Results');
else
    warning('No Stage-1 result files found in %s', Path_Output);
end


%% Plot setting

set(groot, 'defaultAxesFontName', 'Times New Roman');
set(groot, 'defaultAxesFontSize', 12);
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultLineMarkerFaceColor','auto');


%% Plot validation_loss

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
xlabel('$\alpha$');
ylabel('Validation Loss');
legend('show', 'Location', 'best', 'Box', 'off');
grid on;

set(gca, 'LooseInset', get(gca, 'TightInset'));

alpha_min = min(Results.alpha);
alpha_max = max(Results.alpha);
xticks(linspace(alpha_min, alpha_max, 16));

% Output
out_png = sprintf('plot_validation_loss_beta_1.png');
saveas(gcf, fullfile(Path_Output, out_png));


%% Plot log-likelihood (Stage 1)

figure('Position',[100 100 700 450]);
hold on;

L_values = unique(Results.L);

for i = 1:length(L_values)
    L_i = L_values(i);
    idx = Results.L == L_i;

    [alpha_sorted, order] = sort(Results.alpha(idx));
    loglik_sorted = Results.logLikStage1(idx);
    loglik_sorted = loglik_sorted(order);

    plot(alpha_sorted, loglik_sorted, '--', 'LineWidth', 1.8, ...
        'DisplayName', sprintf('L = %d', L_i));
end

hold off;
xlabel('$\alpha$');
ylabel('Stage-1 Log-Likelihood');
legend('show', 'Location', 'best', 'Box', 'off');
grid on;

set(gca, 'LooseInset', get(gca, 'TightInset'));

alpha_min = min(Results.alpha);
alpha_max = max(Results.alpha);
xticks(linspace(alpha_min, alpha_max, 16));

% Output
out_png = 'plot_logLikStage1_beta_1.png';
saveas(gcf, fullfile(Path_Output, out_png));


%% Plot M curve

% Find the smallest validation_loss
best_rows = cell(3,1);
for L_i = 1:3
    rows = Results(Results.L==L_i & Results.beta==1, :);
    [~, ix] = min(rows.validation_loss);
    best_rows{L_i} = rows(ix, :);
end

% sample split
Realized_Return_front = Realized_Return(1:T1, :);
Realized_Return_back = Realized_Return(idx_valid, :);
Realized_Return_all  = Realized_Return;

Risk_Free_Rate_front  = Risk_Free_Rate(1:T1);
Risk_Free_Rate_back  = Risk_Free_Rate(idx_valid);
Risk_Free_Rate_all  = Risk_Free_Rate(:);

dates       = Realized_Return.date;
dates_front = dates(1:T1);
dates_back  = dates(idx_valid);
dates_all   = dates(1:T);

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

% Add paths
Path_Code_08 = fullfile(Path_MainFolder, 'Code', ...
    '08  MLE -  Exponential Polynomial with Distortion');
addpath(Path_Code_08);

function [R_axis, M_bar] = compute_M_curve(gamma_hat, L, alpha, beta, ...
    R_vec, Rf_vec, dates_vec, Smooth_AllR, Smooth_AllR_RND)

    if istable(R_vec),     R_vec     = table2array(R_vec);     end
    if istable(Rf_vec),    Rf_vec    = table2array(Rf_vec);    end
    if istable(dates_vec), dates_vec = table2array(dates_vec); end

    [~, delta_vec, M_vec] = log_likelihood_function( ...
        gamma_hat, R_vec, Rf_vec, L, ...
        Smooth_AllR, Smooth_AllR_RND, dates_vec, true, alpha, beta); %#ok<ASGLU>

    R_axis = Smooth_AllR.(num2str(dates_vec(1)));
    R_base = R_axis(:)';
    N      = numel(R_base);

    date_fields = arrayfun(@(d) num2str(d), dates_vec, 'UniformOutput', false);
    T = numel(dates_vec);
    M_interp = NaN(T, N);
    for t = 1:T
        R_t = Smooth_AllR.(date_fields{t});
        M_t = M_vec(t, :);
        M_interp(t, :) = interp1(R_t, M_t, R_base, 'pchip', NaN);
    end
    M_bar  = mean(M_interp, 1, 'omitnan');    
end


% plot
figure('Position',[50 80 1100 380]);
layout = tiledlayout(1, 3, 'TileSpacing', 'Compact', 'Padding', 'None');

for sp = 1:3
    nexttile;
    hold on;
    for L_i = 1:3
        row = best_rows{L_i};

        Sfile = fullfile(Path_Output, row.file);
        data  = load(Sfile, 'gamma_hat', 'L', 'alpha', 'beta');
        fprintf('\n[Sample %d]  L = %d, alpha = %.1f, beta = %.1f, gamma_hat = %s\n', ...
            sp, row.L, row.alpha, row.beta, mat2str(data.gamma_hat, 4));

        [R_axis, M_bar] = compute_M_curve( ...
            data.gamma_hat, row.L, row.alpha, row.beta, ...
            samples{sp}.R, samples{sp}.Rf, samples{sp}.dates, ...
            Smooth_AllR, Smooth_AllR_RND);

        plot(R_axis, M_bar, 'LineWidth', 1.8, ...
            'DisplayName', sprintf('$L=%d$', row.L));
    end
    hold off;
    title(samples{sp}.name);
    xlabel('$R$'); ylabel('$E(M)$');
    legend('show','Location','northeast','Box','off');
    grid on;
    xlim([0.8 1.2]);
    ylim([0.4 2.3]);
    set(gca,'LooseInset',get(gca,'TightInset'));
end

% Output
out_png = 'plot_M_curve_beta_1.png';
saveas(gcf, fullfile(Path_Output, out_png));

