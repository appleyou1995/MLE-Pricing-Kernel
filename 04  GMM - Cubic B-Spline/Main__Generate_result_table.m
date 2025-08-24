clear; clc;

Path_MainFolder = 'D:\Google\我的雲端硬碟\學術｜研究與論文\論文著作\MLE Pricing Kernel';
Path_Data_04 = fullfile(Path_MainFolder, 'Code', '04  Output');


%% Read estimation reult

b_list = [4,6,8];
Results = struct();

for b = b_list
    fname = sprintf('GMM_theta_b%d.mat', b);
    S = load(fullfile(Path_Data_04, fname));
    Results.(sprintf('b%d', b)) = S;
end

clear b S


%% Build table

max_b = max(b_list);
param_matrix = NaN(max_b+2, numel(b_list));

row_labels = cell(max_b+2,1);
row_labels{1} = 'Jval';

% Fill row labels
for i = 0:max_b
    row_labels{i+2} = sprintf('\\theta_{%d}', i);
end

% Fill value
for j = 1:numel(b_list)
    b = b_list(j);
    field = sprintf('b%d', b);

    % Log-likelihood / Jval
    param_matrix(1,j) = Results.(field).Jval;

    % Theta values
    theta_hat = Results.(field).theta_hat;
    param_matrix(2:(b+2), j) = theta_hat(:);
end

T = array2table(param_matrix, 'VariableNames', ...
    strcat('b', string(b_list)), 'RowNames', row_labels);

disp(T)