clear; clc;
Path_MainFolder = 'D:\Google\我的雲端硬碟\學術｜研究與論文\論文著作\MLE Pricing Kernel';
Path_Output = fullfile(Path_MainFolder, 'Code', '11  Output');


%% Summarize MLE Results (Loglik & BIC)

folder = Path_Output;
files  = dir(fullfile(folder, 'MLE_gamma_L_*.mat'));

S = containers.Map('KeyType','char','ValueType','any');

init_row = @() struct( ...
    'alpha', NaN, 'beta', NaN, ...
    'L1_loglik', NaN, 'L1_BIC', NaN, ...
    'L2_loglik', NaN, 'L2_BIC', NaN, ...
    'L3_loglik', NaN, 'L3_BIC', NaN);

fprintf('Processing %d files...\n', numel(files));

for k = 1:numel(files)
    try
        data = load(fullfile(files(k).folder, files(k).name), ...
                    'alpha','beta','L','log_lik','BIC');
    catch
        warning('Unable to load file: %s', files(k).name);
        continue;
    end
    
    a = data.alpha;
    b = data.beta;
    L = data.L;
    
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
            row.L1_loglik = data.log_lik;
            row.L1_BIC    = data.BIC;
        case 2
            row.L2_loglik = data.log_lik;
            row.L2_BIC    = data.BIC;
        case 3
            row.L3_loglik = data.log_lik;
            row.L3_BIC    = data.BIC;
        otherwise
            warning('Warning: L=%d in file %s', L, files(k).name);
    end
    
    S(key) = row;
end

rows = values(S);
rows = [rows{:}];
T = struct2table(rows);

T.alpha = round(T.alpha, 2);
T.beta  = round(T.beta,  2);

vars = T.Properties.VariableNames;
for i = 1:numel(vars)
    if isnumeric(T.(vars{i})) && ~strcmp(vars{i}, 'alpha') && ~strcmp(vars{i}, 'beta')
        T.(vars{i}) = round(T.(vars{i}), 4);
    end
end

T = sortrows(T, {'alpha','beta'});


%% Find Best Models (Max Loglik & Min BIC)

% === 1. Max Log-Likelihood ===
[max_LL1, idx_L1] = max(T.L1_loglik, [], 'omitnan');
[max_LL2, idx_L2] = max(T.L2_loglik, [], 'omitnan');
[max_LL3, idx_L3] = max(T.L3_loglik, [], 'omitnan');

% === 2. Min BIC ===
[min_BIC1, idx_B1] = min(T.L1_BIC, [], 'omitnan');
[min_BIC2, idx_B2] = min(T.L2_BIC, [], 'omitnan');
[min_BIC3, idx_B3] = min(T.L3_BIC, [], 'omitnan');

clc;
fprintf('\n======================================================\n');
fprintf('                 MODEL SELECTION SUMMARY\n');
fprintf('======================================================\n');

fprintf('\n--- Max Log-Likelihood (Higher is Better) ---\n');
if ~isnan(max_LL1)
    fprintf('L = 1: Loglik = %8.4f | (alpha=%.2f, beta=%.2f)\n', max_LL1, T.alpha(idx_L1), T.beta(idx_L1));
else
    fprintf('L = 1: No Data\n');
end
if ~isnan(max_LL2)
    fprintf('L = 2: Loglik = %8.4f | (alpha=%.2f, beta=%.2f)\n', max_LL2, T.alpha(idx_L2), T.beta(idx_L2));
else
    fprintf('L = 2: No Data\n');
end
if ~isnan(max_LL3)
    fprintf('L = 3: Loglik = %8.4f | (alpha=%.2f, beta=%.2f)\n', max_LL3, T.alpha(idx_L3), T.beta(idx_L3));
else
    fprintf('L = 3: No Data\n');
end

fprintf('\n--- Min BIC (Lower is Better) ---\n');
if ~isnan(min_BIC1)
    fprintf('L = 1: BIC    = %8.4f | (alpha=%.2f, beta=%.2f)\n', min_BIC1, T.alpha(idx_B1), T.beta(idx_B1));
else
    fprintf('L = 1: No Data\n');
end
if ~isnan(min_BIC2)
    fprintf('L = 2: BIC    = %8.4f | (alpha=%.2f, beta=%.2f)\n', min_BIC2, T.alpha(idx_B2), T.beta(idx_B2));
else
    fprintf('L = 2: No Data\n');
end
if ~isnan(min_BIC3)
    fprintf('L = 3: BIC    = %8.4f | (alpha=%.2f, beta=%.2f)\n', min_BIC3, T.alpha(idx_B3), T.beta(idx_B3));
else
    fprintf('L = 3: No Data\n');
end
fprintf('======================================================\n');


%% Output csv

out_csv = fullfile(folder, 'MLE_estimation_summary.csv');
writetable(T, out_csv);
fprintf('\nSummary file saved to:\n%s\n', out_csv);