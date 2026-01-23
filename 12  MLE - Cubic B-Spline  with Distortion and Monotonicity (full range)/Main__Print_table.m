clear; clc;

Path_MainFolder = 'D:\Google\我的雲端硬碟\學術｜研究與論文\論文著作\MLE Pricing Kernel';
Path_Output = fullfile(Path_MainFolder, 'Code', '12  Output');


%% Add MLE log-likelihood

function x = pick(v,i)
    if numel(v) >= i, x = v(i); else, x = NaN; end
end

folder = Path_Output;
files  = dir(fullfile(folder, 'MLE_BSpline_b_*.mat'));

S = containers.Map('KeyType','char','ValueType','any');

% 初始化欄位設定
% b=4 -> 5 params; 
% b=6 -> 7 params; 
% b=8 -> 9 params
init_row = @() struct( ...
    'alpha', NaN, 'beta', NaN, ...
    'b4_loglik', NaN, ...
    'b4_theta1', NaN, 'b4_theta2', NaN, 'b4_theta3', NaN, 'b4_theta4', NaN, 'b4_theta5', NaN, ...
    'b6_loglik', NaN, ...
    'b6_theta1', NaN, 'b6_theta2', NaN, 'b6_theta3', NaN, 'b6_theta4', NaN, 'b6_theta5', NaN, 'b6_theta6', NaN, 'b6_theta7', NaN, ...
    'b8_loglik', NaN, ...
    'b8_theta1', NaN, 'b8_theta2', NaN, 'b8_theta3', NaN, 'b8_theta4', NaN, 'b8_theta5', NaN, 'b8_theta6', NaN, 'b8_theta7', NaN, 'b8_theta8', NaN, 'b8_theta9', NaN);

for k = 1:numel(files)
    % 讀取 b_val 與 theta_hat
    try
        data = load(fullfile(files(k).folder, files(k).name), ...
                    'alpha','beta','b_val','theta_hat','log_lik');
    catch
        warning('無法讀取檔案: %s', files(k).name);
        continue;
    end
    
    a = data.alpha;
    b_dist = data.beta;
    b_spline = data.b_val;
    t = data.theta_hat(:)';
    
    key = sprintf('a%.2f_b%.2f', a, b_dist);
    
    if S.isKey(key)
        row = S(key);
    else
        row = init_row();
        row.alpha = a;
        row.beta  = b_dist;
    end
    
    % 根據 b_spline 填入對應欄位
    switch b_spline
        case 4
            row.b4_loglik = data.log_lik;
            for i = 1:5, row.(['b4_theta' num2str(i)]) = pick(t, i); end
            
        case 6
            row.b6_loglik = data.log_lik;
            for i = 1:7, row.(['b6_theta' num2str(i)]) = pick(t, i); end
            
        case 8
            row.b8_loglik = data.log_lik;
            for i = 1:9, row.(['b8_theta' num2str(i)]) = pick(t, i); end
            
        otherwise
            warning('Warning: Unexpected b=%d in file %s', b_spline, files(k).name);
    end
    
    S(key) = row;
end

clear k data a b_dist b_spline t key row files

rows = values(S);
if isempty(rows)
    error('沒有找到任何符合的資料，請檢查路徑或檔案是否已產生。');
end
rows = [rows{:}];
T = struct2table(rows);

% 格式化數值
T.alpha = round(T.alpha, 2);
T.beta  = round(T.beta,  2);

vars = T.Properties.VariableNames;
for i = 1:numel(vars)
    col = T.(vars{i});
    if isnumeric(col)
        T.(vars{i}) = round(col, 4);
    end
end

% 排序
T = sortrows(T, {'alpha','beta'});

% 修改 5: 找出各個 b 的最佳 LogLik
[best4, idx4] = max(T.b4_loglik, [], 'omitnan');
[best6, idx6] = max(T.b6_loglik, [], 'omitnan');
[best8, idx8] = max(T.b8_loglik, [], 'omitnan');

% 防止空資料報錯
if isempty(best4), best4 = -Inf; idx4 = 1; end
if isempty(best6), best6 = -Inf; idx6 = 1; end
if isempty(best8), best8 = -Inf; idx8 = 1; end

best4_alpha = T.alpha(idx4); best4_beta = T.beta(idx4);
best6_alpha = T.alpha(idx6); best6_beta = T.beta(idx6);
best8_alpha = T.alpha(idx8); best8_beta = T.beta(idx8);

clc
fprintf('\n=== Max loglik by b (Spline Knots) ===\n');
fprintf('b = 4: alpha = %.2f, beta = %.2f, loglik = %.4f\n', best4_alpha, best4_beta, best4);
fprintf('b = 6: alpha = %.2f, beta = %.2f, loglik = %.4f\n', best6_alpha, best6_beta, best6);
fprintf('b = 8: alpha = %.2f, beta = %.2f, loglik = %.4f\n', best8_alpha, best8_beta, best8);


%% Output csv

out_csv = fullfile(folder, 'MLE_BSpline_estimation_summary.csv');
writetable(T, out_csv);
fprintf('Summary saved to: %s\n', out_csv);