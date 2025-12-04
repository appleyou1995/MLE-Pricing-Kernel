clear; clc;

%% Path setting

Path_MainFolder = 'D:\Google\我的雲端硬碟\學術｜研究與論文\論文著作';

Path_TTM = fullfile(Path_MainFolder, 'CDI Method',         'Code', '01  輸出資料');
Path_MLE = fullfile(Path_MainFolder, 'MLE Pricing Kernel', 'Code', '10  Output');

Path_Output = Path_MLE;


%% List setting

TTM_list = [30, 60, 90, 180];
L_fixed  = 1;


%% Part 1: kappa table

Kappa_All = table();

for i = 1:numel(TTM_list)

    TTM = TTM_list(i);

    % 1. Load TTM_XX.csv
    File_TTM = fullfile(Path_TTM, sprintf('TTM_%d.csv', TTM));
    if ~isfile(File_TTM)
        warning('TTM file not found: %s', File_TTM);
        continue;
    end

    T_tmp = readtable(File_TTM);

    % 2. Load kappa_vec
    File_MLE = fullfile(Path_MLE, sprintf('MLE_gamma_TTM_%d_L_%d.mat', TTM, L_fixed));
    if ~isfile(File_MLE)
        warning('MLE file not found: %s', File_MLE);
        continue;
    end

    S = load(File_MLE);

    if ~isfield(S, 'kappa_vec')
        error('File %s has no kappa_vec.', File_MLE);
    end

    kappa_vec = S.kappa_vec(:);    % T×1

    % 3. Length checking
    n_csv   = height(T_tmp);
    n_kappa = numel(kappa_vec);
    
    if n_csv ~= n_kappa
        warning('Length mismatch for TTM = %d: CSV rows = %d, kappa_vec length = %d → trimming to match.', ...
            TTM, n_csv, n_kappa);
    
        if n_csv > n_kappa
            T_tmp = T_tmp(1:n_kappa, :);
        else
            kappa_vec = kappa_vec(1:n_csv);
        end
    end

    % 4. Rename column & add kappa
    %    date  → date_TTM_XX
    %    exdate → exdate_TTM_XX
    new_date_var   = sprintf('date_TTM_%d',   TTM);
    new_exdate_var = sprintf('exdate_TTM_%d', TTM);
    new_kappa_var  = sprintf('kappa_TTM_%d',  TTM);

    if any(strcmp('date',   T_tmp.Properties.VariableNames)) && ...
       any(strcmp('exdate', T_tmp.Properties.VariableNames))

        T_tmp = renamevars(T_tmp, {'date','exdate'}, ...
                                   {new_date_var, new_exdate_var});
    else
        error('TTM_%d.csv does not contain expected variables "date" and "exdate".', TTM);
    end

    T_tmp.(new_kappa_var) = kappa_vec;

    % 5. Merge
    if isempty(Kappa_All)
        Kappa_All = T_tmp;
    else
        Kappa_All = [Kappa_All, T_tmp];  %#ok<AGROW>
    end

end

% 6. Output CSV
File_KappaCSV = fullfile(Path_Output, 'Kappa_By_TTM.csv');
writetable(Kappa_All, File_KappaCSV);
fprintf('Kappa table saved to: %s\n', File_KappaCSV);


%% Part 2: gamma table

TTM_col   = [];
gamma_col = [];

for i = 1:numel(TTM_list)

    TTM = TTM_list(i);

    File_MLE = fullfile(Path_MLE, sprintf('MLE_gamma_TTM_%d_L_%d.mat', TTM, L_fixed));
    if ~isfile(File_MLE)
        warning('MLE file not found: %s', File_MLE);
        continue;
    end

    S = load(File_MLE);

    if ~isfield(S, 'gamma_hat')
        error('File %s has no gamma_hat.', File_MLE);
    end

    g = S.gamma_hat(:);

    if numel(g) ~= 1
        warning('gamma_hat for TTM = %d has length %d; using first element.', TTM, numel(g));
        g = g(1);
    end

    TTM_col   = [TTM_col;   TTM];  %#ok<AGROW>
    gamma_col = [gamma_col; g];    %#ok<AGROW>

end

Gamma_Table = table(TTM_col, gamma_col, ...
    'VariableNames', {'TTM', 'gamma'});

File_GammaCSV = fullfile(Path_Output, 'Gamma_By_TTM.csv');
writetable(Gamma_Table, File_GammaCSV);
fprintf('Gamma table saved to: %s\n', File_GammaCSV);


%% Plot kappa

set(groot, 'defaultAxesFontName', 'Times New Roman');
set(groot, 'defaultAxesFontSize', 12);
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');

t30  = datetime(Kappa_All.date_TTM_30,  'ConvertFrom','yyyymmdd');
t60  = datetime(Kappa_All.date_TTM_60,  'ConvertFrom','yyyymmdd');
t90  = datetime(Kappa_All.date_TTM_90,  'ConvertFrom','yyyymmdd');
t180 = datetime(Kappa_All.date_TTM_180, 'ConvertFrom','yyyymmdd');

k30  = Kappa_All.kappa_TTM_30;
k60  = Kappa_All.kappa_TTM_60;
k90  = Kappa_All.kappa_TTM_90;
k180 = Kappa_All.kappa_TTM_180;

fig = figure('Position',[200 200 850 400]);
tiledlayout(1, 1, 'TileSpacing', 'Compact', 'Padding', 'None');
nexttile;
hold on; grid on;

plot(t30,  k30,  'LineWidth',1.3);
plot(t60,  k60,  'LineWidth',1.3);
plot(t90,  k90,  'LineWidth',1.3);
plot(t180, k180, 'LineWidth',1.3);

ylabel('$\kappa_t$',Rotation=0);

legend({'TTM = 30', 'TTM = 60', 'TTM = 90', 'TTM = 180'}, ...
       'Location','best','Box','off');

hold off;
saveas(fig, fullfile(Path_Output, 'Figure_kappa.png'));