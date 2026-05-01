clear; clc;

Path_Main = 'D:\Google\我的雲端硬碟\學術｜研究與論文\論文著作\MLE Pricing Kernel\Code';
Path_Output = fullfile(Path_Main, '16  Output - Summary Table');


%% Master Switch

Target_TTM = 30;
% [30, 60, 90, 180]

Extension_Mode = 'Full';
% ['Full', 'Non_Recession', 'High_Vol', 'Low_Vol']

fprintf('Current Mode: %s (TTM = %d)\n', Extension_Mode, Target_TTM);


%% Extension_Mode

switch Extension_Mode
    case 'Full'
        Path_Result = fullfile(Path_Main, '16  Output (TTM Comparison)');
    case 'Non_Recession'
        Path_Result = fullfile(Path_Main, '16  Output (Non-Recession Periods)');
    case 'High_Vol'
        Path_Result = fullfile(Path_Main, '16  Output (High vs. Low Volatility)');
    case 'Low_Vol'
        Path_Result = fullfile(Path_Main, '16  Output (High vs. Low Volatility)');
end


%% Summarize MLE Results

folder = Path_Result;
file_pattern = sprintf('MLE_BSpline_TTM_%d_%s_*.mat', Target_TTM, Extension_Mode);
files = dir(fullfile(folder, file_pattern));

if isempty(files)
    error('在資料夾 %s 中找不到符合條件的檔案 (%s)，請檢查是否已執行估計。', folder, file_pattern);
end

S = containers.Map('KeyType','char','ValueType','any');
b_range = 6; 
init_struct = struct('alpha', NaN, 'beta', NaN);
for b = b_range
    init_struct.(sprintf('b%d_loglik', b)) = NaN;
end
init_row = @() init_struct;

fprintf('Processing %d files in %s...\n', numel(files), Extension_Mode);

for k = 1:numel(files)
    try
        data = load(fullfile(files(k).folder, files(k).name), ...
                    'alpha','beta','b_val','log_lik');
    catch
        warning('無法讀取檔案: %s', files(k).name);
        continue;
    end
    
    alpha_dist = data.alpha;
    beta_dist = data.beta;
    b_spline = data.b_val;
    
    key = sprintf('a%.2f_b%.2f', alpha_dist, beta_dist);
    
    if S.isKey(key)
        row = S(key);
    else
        row = init_row();
        row.alpha = alpha_dist;
        row.beta  = beta_dist;
    end
    
    field_ll  = sprintf('b%d_loglik', b_spline);    
    row.(field_ll)  = data.log_lik;
    
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


%% Find Best Models

clc;
fprintf('\n==============================================================\n');
fprintf('     B-SPLINE MODEL SELECTION SUMMARY (%s, TTM=%d)\n', Extension_Mode, Target_TTM);
fprintf('--------------------------------------------------------------\n');

% Max Log-Likelihood
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

fprintf('==============================================================\n');


%% Transpose

target_b = 6;
ll_field = sprintf('b%d_loglik', target_b);
T_matrix_ll = unstack(T(:, {'alpha', 'beta', ll_field}), ...
    ll_field, 'beta', ...
    'VariableNamingRule', 'preserve');
old_names = T_matrix_ll.Properties.VariableNames;
new_names = old_names;
for v = 2:numel(old_names)
    num_str = strrep(old_names{v}, 'x', '');
    num_str = strrep(num_str, '_', '.');

    val = str2double(num_str);
    new_names{v} = sprintf('beta_%.2f', val);
end
T_matrix_ll.Properties.VariableNames = new_names;

alpha_strs = arrayfun(@(x) sprintf('%.2f', x), T_matrix_ll.alpha, 'UniformOutput', false);
T_matrix_ll.alpha = alpha_strs;

out_csv_name = sprintf('Summary_TTM_%d_%s.csv', Target_TTM, Extension_Mode);
out_csv = fullfile(Path_Output, out_csv_name);

csv_headers = T_matrix_ll.Properties.VariableNames;
csv_headers{1} = 'alpha';
T_matrix_ll.Properties.VariableNames = csv_headers;
disp(T_matrix_ll);


%% Output csv

writetable(T_matrix_ll, out_csv);
fprintf('\nSummary file saved to:\n%s\n', out_csv);