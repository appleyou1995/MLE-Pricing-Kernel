clear; clc;
Path_MainFolder = 'D:\Google\我的雲端硬碟\學術｜研究與論文\論文著作\MLE Pricing Kernel';
% Path_Output = fullfile(Path_MainFolder, 'Code', '12  Output');
Path_Output = fullfile(Path_MainFolder, 'Code', '12  Output - Quartic');


%% Summarize MLE Results

folder  = Path_Output;
files   = dir(fullfile(folder, 'MLE_BSpline_b_*.mat'));
S       = containers.Map('KeyType','char','ValueType','any');
b_range = 4:9; 

init_struct = struct('alpha', NaN, 'beta', NaN);
for b = b_range
    init_struct.(sprintf('b%d_loglik', b)) = NaN;
    init_struct.(sprintf('b%d_BIC', b))    = NaN;
end
init_row = @() init_struct;

fprintf('Processing %d files...\n', numel(files));

for k = 1:numel(files)
    try
        data = load(fullfile(files(k).folder, files(k).name), ...
                    'alpha','beta','b_val','log_lik','BIC');
    catch
        warning('無法讀取檔案: %s', files(k).name);
        continue;
    end
    
    a = data.alpha;
    b_dist = data.beta;
    b_spline = data.b_val;
    
    key = sprintf('a%.2f_b%.2f', a, b_dist);
    
    if S.isKey(key)
        row = S(key);
    else
        row = init_row();
        row.alpha = a;
        row.beta  = b_dist;
    end
    
    field_ll  = sprintf('b%d_loglik', b_spline);
    field_bic = sprintf('b%d_BIC', b_spline);
    
    row.(field_ll)  = data.log_lik;
    row.(field_bic) = data.BIC;
    
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

clc;
fprintf('\n======================================================\n');
fprintf('           B-SPLINE MODEL SELECTION SUMMARY\n');
fprintf('======================================================\n');

% === 1. Max Log-Likelihood (Higher is Better) ===
fprintf('\n--- Max Log-Likelihood (Higher is Better) ---\n');
for b = b_range
    col_name = sprintf('b%d_loglik', b);
    if ismember(col_name, T.Properties.VariableNames)
        [max_val, idx] = max(T.(col_name), [], 'omitnan');
        if ~isempty(max_val) && ~isnan(max_val)
            fprintf('b = %d: Loglik = %8.4f | (alpha=%.2f, beta=%.2f)\n', ...
                b, max_val, T.alpha(idx), T.beta(idx));
        else
            fprintf('b = %d: No Data\n', b);
        end
    end
end

% === 2. Min BIC (Lower is Better) ===
fprintf('\n--- Min BIC (Lower is Better) ---\n');
for b = b_range
    col_name = sprintf('b%d_BIC', b);
    if ismember(col_name, T.Properties.VariableNames)
        [min_val, idx] = min(T.(col_name), [], 'omitnan');
        if ~isempty(min_val) && ~isnan(min_val)
            fprintf('b = %d: BIC    = %8.4f | (alpha=%.2f, beta=%.2f)\n', ...
                b, min_val, T.alpha(idx), T.beta(idx));
        else
            fprintf('b = %d: No Data\n', b);
        end
    end
end
fprintf('======================================================\n');


%% Output csv

out_csv = fullfile(folder, 'MLE_BSpline_estimation_summary.csv');
writetable(T, out_csv);
fprintf('\nSummary file saved to:\n%s\n', out_csv);