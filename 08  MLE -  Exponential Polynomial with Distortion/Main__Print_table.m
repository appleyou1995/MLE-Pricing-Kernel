clc

% ----- Case 1 -----
% Path_Output = fullfile(Path_MainFolder, 'Code', '08  Output - half');

% ----- Case 2 -----
Path_Output = fullfile(Path_MainFolder, 'Code', '08  Output - two_thirds');

% ----- Case 3 -----
% Path_Output = fullfile(Path_MainFolder, 'Code', '08  Output - all');


%% 

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

    a = data.alpha; b = data.beta; L = data.L;
    g = data.gamma_hat(:)';
    key = sprintf('a%.6f_b%.6f', a, b);

    if S.isKey(key)
        row = S(key);
    else
        row = init_row();
        row.alpha = a; row.beta = b;
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

rows = values(S);
rows = [rows{:}];
T = struct2table(rows);
T = sortrows(T, {'alpha','beta'});
T.alpha = round(T.alpha, 1);
vars = T.Properties.VariableNames;
for i = 1:numel(vars)
    col = T.(vars{i});
    if isnumeric(col)
        T.(vars{i}) = round(col, 2);
    end
end

% beta = 1
T = T(T.beta == 1, :);

% Output csv
out_csv = fullfile(folder, 'Stage1_MLE_estimation_summary.csv');
writetable(T, out_csv);