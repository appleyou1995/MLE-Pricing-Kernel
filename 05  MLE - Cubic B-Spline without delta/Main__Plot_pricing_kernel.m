clear; clc;

Path_MainFolder = 'D:\Google\我的雲端硬碟\學術｜研究與論文\論文著作\MLE Pricing Kernel';
Path_Output = fullfile(Path_MainFolder, 'Code', '05  Output');


%% Load the data

Path_Data = 'D:\Google\我的雲端硬碟\學術｜研究與論文\論文著作\CDI Method';

Target_TTM = 30;

% Load Q-measure PDF tables: R axis and corresponding f^*_t(R)
Path_Data_02 = fullfile(Path_Data, 'Code', '02  輸出資料');
Smooth_AllR = [];

years_to_merge = 1996:2021;
for year = years_to_merge
    input_filename = fullfile(Path_Data_02, sprintf('TTM_%d_RND_Tables_%d.mat', Target_TTM, year));
    if exist(input_filename, 'file')
        data = load(input_filename);
        Smooth_AllR = [Smooth_AllR, data.Table_Smooth_AllR];               % R_grid for interpolation
    else
        warning('File %s does not exist.', input_filename);
    end
end

clear input_filename year Path_Data Path_Data_02 data years_to_merge


%% Load estimation result

mat_files = dir(fullfile(Path_Output, 'MLE_theta_b*.mat'));

for k = 1:length(mat_files)
    file_path = fullfile(Path_Output, mat_files(k).name);
    load(file_path, 'M_vec');
    b_value = regexp(mat_files(k).name, '(?<=_b)(\d+)', 'match', 'once');
    var_name = ['M_vec_' b_value];
    assignin('base', var_name, M_vec);
end

clear file_path b_value var_name k M_vec mat_files


%% Average M across all periods and Plot

b_list = [4, 6, 8];
M_all = {M_vec_4, M_vec_6, M_vec_8};

date_fields = Smooth_AllR.Properties.VariableNames;
R_base = Smooth_AllR.(date_fields{1});

T = size(date_fields, 2);
N = numel(R_base);

figure;
set(gcf, 'Position', [100, 100, 1200, 400]);
layout = tiledlayout(1, 3, 'TileSpacing','Compact','Padding','None');

for j = 1:length(b_list)
    b = b_list(j);
    M_vec = M_all{j};

    M_interp = NaN(T, N);
    for t = 1:T
        R_t = Smooth_AllR.(date_fields{t});
        M_t = M_vec(t,:);
        M_interp(t,:) = interp1(R_t, M_t, R_base, 'linear', NaN);
    end
    
    M_avg = mean(M_interp, 1, 'omitnan');
    
    nexttile;
    plot(R_base, M_avg, 'LineWidth', 1.5);
    set(gca, 'FontName','Times New Roman');
    grid on;
    xlim([0.8 1.2]);
    %ylim([0.8 2.5]);

    xlabel('Gross Return', 'FontName','Times New Roman');
    if j == 1
        ylabel('$\mathrm{E}(M)$', 'Interpreter','latex');
    end

    title(sprintf('$b=%d$', b), ...
          'FontName','Times New Roman', 'Interpreter','latex');
end

filename = 'Average Pricing Kernel.png';
saveas(gcf, fullfile(Path_Output, filename));


%% Average M across all periods

b_list   = [4, 6, 8];
M_all    = {M_vec_4, M_vec_6, M_vec_8};

date_fields = Smooth_AllR.Properties.VariableNames;
R_base = Smooth_AllR.(date_fields{1});
T = numel(date_fields);
N = numel(R_base);

M_avg_all = NaN(numel(b_list), N);

for j = 1:numel(b_list)
    M_vec = M_all{j};

    M_interp = NaN(T, N);
    for t = 1:T
        R_t = Smooth_AllR.(date_fields{t});
        M_t = M_vec(t,:);
        M_interp(t,:) = interp1(R_t, M_t, R_base, 'linear', NaN);
    end

    M_avg_all(j,:) = mean(M_interp, 1, 'omitnan');
end


%% Plot setting

set(groot, 'defaultAxesFontName', 'Times New Roman');
set(groot, 'defaultAxesFontSize', 12);
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');


%% Plot Average Pricing Kernel

figure; set(gcf, 'Position', [120, 120, 680, 480]);

hold on;
for j = 1:numel(b_list)
    plot(R_base, M_avg_all(j,:), 'LineWidth', 1.5, ...
         'DisplayName', sprintf('$b = %d$', b_list(j)));
end
hold off;

grid on;
box on;
xlim([0.8 1.2]);
%ylim([0.8 2.5]);

xlabel('Gross Return $R$');
ylabel('$\mathrm{E}(M)$');
title('Average Pricing Kernel (MLE)');
legend('Location','southwest', 'Box','off');

saveas(gcf, fullfile(Path_Output, 'Average_Pricing_Kernel_MLE.png'));