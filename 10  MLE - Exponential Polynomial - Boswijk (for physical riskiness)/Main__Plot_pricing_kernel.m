clear; clc;

Path_MainFolder = 'D:\Google\我的雲端硬碟\學術｜研究與論文\論文著作\MLE Pricing Kernel';
Path_Data = 'D:\Google\我的雲端硬碟\學術｜研究與論文\論文著作\CDI Method';
Path_Output = fullfile(Path_MainFolder, 'Code', '10  Output');


%% Plot setting

set(groot, 'defaultAxesFontName', 'Times New Roman');
set(groot, 'defaultAxesFontSize', 12);
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');


%% Grid setting

x_min = 0.8;
x_max = 1.2;

TTM_list = [30, 60, 90, 180];
R_axis = (x_min : 0.001 : x_max)';
N = numel(R_axis);

M_avg_allTTM = NaN(numel(TTM_list), N);


%% Load estimation result and construct M curve

for i = 1:numel(TTM_list)

    TTM = TTM_list(i);
    FileName = fullfile(Path_Output, sprintf('MLE_gamma_TTM_%d_L_%d.mat', TTM, 1));

    if ~isfile(FileName)
        warning('File not found: %s', FileName);
        continue;
    end

    S = load(FileName);

    gamma_hat = S.gamma_hat(:);                                            % L×1
    kappa_vec = S.kappa_vec(:);                                            % T×1
    L         = numel(gamma_hat);
    T         = numel(kappa_vec);

    % ---- Step 1 ----
    logR = log(R_axis);                                                    % N×1
    poly_sum = zeros(N,1);
    for l = 1:L
        poly_sum = poly_sum + gamma_hat(l) .* (logR.^l);
    end
    % poly_sum: N×1

    % ---- Step 2 ----
    % logM_tR = κ_t - poly_sum(R)
    logM_all = kappa_vec * ones(1, N) - ones(T,1) * poly_sum.';            % T×N
    M_all    = exp(logM_all);                                              % T×N

    % ---- Step 3 ----
    M_avg = mean(M_all, 1);                                                % 1×N
    M_avg_allTTM(i, :) = M_avg;

end


%% Plot

fig = figure('Position',[200 200 550 450]);
tiledlayout(1, 1, 'TileSpacing', 'Compact', 'Padding', 'None');
nexttile;
hold on; grid on;

colors = lines(numel(TTM_list));

for i = 1:numel(TTM_list)
    if all(isnan(M_avg_allTTM(i,:)))
        continue;
    end
    plot(R_axis, M_avg_allTTM(i, :), 'LineWidth', 1.5, ...
        'Color', colors(i,:));
end

xlim([x_min x_max]);
ylim([0.65 1.65]);
xlabel('Gross return $R_{t+1}$');
ylabel('Average pricing kernel $\bar{M}(R)$');
legend({'TTM = 30', 'TTM = 60', 'TTM = 90', 'TTM = 180'}, ...
       'Location','best','Box','off');

hold off;

saveas(fig, fullfile(Path_Output, 'Figure_Avg_Pricing_Kernel.png'));