clear; clc;

Path_MainFolder = 'D:\Google\我的雲端硬碟\學術｜研究與論文\論文著作\MLE Pricing Kernel';
Path_Output = fullfile(Path_MainFolder, 'Code', '09  Output');


%% Add MLE log-likelihood

function x = pick(v,i)
    if numel(v) >= i, x = v(i); else, x = NaN; end
end

folder = Path_Output;
files  = dir(fullfile(folder, 'MLE_gamma_L_*.mat'));

S = containers.Map('KeyType','char','ValueType','any');

init_row = @() struct( ...
    'alpha', NaN, 'beta', NaN, ...
    'L1_gamma1', NaN, 'L1_loglik', NaN, ...
    'L2_gamma1', NaN, 'L2_gamma2', NaN, 'L2_loglik', NaN, ...
    'L3_gamma1', NaN, 'L3_gamma2', NaN, 'L3_gamma3', NaN, 'L3_loglik', NaN);


for k = 1:numel(files)
    data = load(fullfile(files(k).folder, files(k).name), ...
                'alpha','beta','L','gamma_hat','log_lik');

    a = data.alpha;
    b = data.beta;
    L = data.L;
    g = data.gamma_hat(:)';
    key = sprintf('a%.2f_b%.2f', a, b);

    if S.isKey(key)
        row = S(key);
    else
        row = init_row();
        row.alpha = a;
        row.beta  = b;
    end

    switch L
        case 1
            row.L1_gamma1 = pick(g,1);
            row.L1_loglik = data.log_lik;
        case 2
            row.L2_gamma1 = pick(g,1);
            row.L2_gamma2 = pick(g,2);
            row.L2_loglik = data.log_lik;
        case 3
            row.L3_gamma1 = pick(g,1);
            row.L3_gamma2 = pick(g,2);
            row.L3_gamma3 = pick(g,3);
            row.L3_loglik = data.log_lik;
        otherwise
            warning('Warning: L=%d（%s）', L, files(k).name);
    end

    S(key) = row;
end

clear k data a b L key row files

rows = values(S);
rows = [rows{:}];
T = struct2table(rows);
T = sortrows(T, {'alpha','beta'});
T.alpha = round(T.alpha, 2);
T.beta  = round(T.beta,  2);
vars = T.Properties.VariableNames;
for i = 1:numel(vars)
    col = T.(vars{i});
    if isnumeric(col)
        T.(vars{i}) = round(col, 2);
    end
end

% beta = 1
% T = T(T.beta == 1, :);

[best1, idx1] = max(T.L1_loglik, [], 'omitnan');
[best2, idx2] = max(T.L2_loglik, [], 'omitnan');
[best3, idx3] = max(T.L3_loglik, [], 'omitnan');

best1_alpha = T.alpha(idx1); best1_beta = T.beta(idx1);
best2_alpha = T.alpha(idx2); best2_beta = T.beta(idx2);
best3_alpha = T.alpha(idx3); best3_beta = T.beta(idx3);

clc
fprintf('\n=== Max loglik by L ===\n');
fprintf('L = 1: alpha = %.2f, beta = %.2f, loglik = %.2f\n', best1_alpha, best1_beta, best1);
fprintf('L = 2: alpha = %.2f, beta = %.2f, loglik = %.2f\n', best2_alpha, best2_beta, best2);
fprintf('L = 3: alpha = %.2f, beta = %.2f, loglik = %.2f\n', best3_alpha, best3_beta, best3);


%% Add GMM validation loss

val_path = fullfile(folder, 'GMM_validation_results.mat');
S2 = load(val_path);
Results = S2.Results;
Results = sortrows(Results, {'L','alpha','beta'});

p = 2;
T.alpha       = round(T.alpha, p);
T.beta        = round(T.beta,  p);
Results.alpha = round(Results.alpha, p);
Results.beta  = round(Results.beta,  p);

T.L1_validation_loss = NaN(height(T), 1);
T.L2_validation_loss = NaN(height(T), 1);
T.L3_validation_loss = NaN(height(T), 1);

fixed_cols = T.Properties.VariableNames(1:2);
other_cols = sort(T.Properties.VariableNames(3:end));
T = T(:, [fixed_cols, other_cols]);

for i = 1:height(Results)
    L_i     = Results.L(i);
    a_i     = Results.alpha(i);
    b_i     = Results.beta(i);
    v_loss  = Results.validation_loss(i);

    idx = (T.alpha == a_i) & (T.beta == b_i);

    switch L_i
        case 1
            T.L1_validation_loss(idx) = v_loss;
        case 2
            T.L2_validation_loss(idx) = v_loss;
        case 3
            T.L3_validation_loss(idx) = v_loss;
    end
end

T.L1_validation_loss = round(T.L1_validation_loss, 8);
T.L2_validation_loss = round(T.L2_validation_loss, 8);
T.L3_validation_loss = round(T.L3_validation_loss, 8);

[best1_val, idx1_val] = min(T.L1_validation_loss, [], 'omitnan');
[best2_val, idx2_val] = min(T.L2_validation_loss, [], 'omitnan');
[best3_val, idx3_val] = min(T.L3_validation_loss, [], 'omitnan');

best1_alpha_val = T.alpha(idx1_val); best1_beta_val = T.beta(idx1_val);
best2_alpha_val = T.alpha(idx2_val); best2_beta_val = T.beta(idx2_val);
best3_alpha_val = T.alpha(idx3_val); best3_beta_val = T.beta(idx3_val);

fprintf('\n=== Min validation loss by L ===\n');
fprintf('L = 1: alpha = %.2f, beta = %.2f, validation_loss = %.8f\n', ...
    best1_alpha_val, best1_beta_val, best1_val);
fprintf('L = 2: alpha = %.2f, beta = %.2f, validation_loss = %.8f\n', ...
    best2_alpha_val, best2_beta_val, best2_val);
fprintf('L = 3: alpha = %.2f, beta = %.2f, validation_loss = %.8f\n', ...
    best3_alpha_val, best3_beta_val, best3_val);


%% Output csv

fmt = @(x) arrayfun(@(y) sprintf('%.2e', y), x, 'UniformOutput', false);

T.L1_validation_loss = fmt(T.L1_validation_loss);
T.L2_validation_loss = fmt(T.L2_validation_loss);
T.L3_validation_loss = fmt(T.L3_validation_loss);

out_csv = fullfile(folder, 'MLE_estimation_summary.csv');
writetable(T, out_csv);