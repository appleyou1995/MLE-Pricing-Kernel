clear; clc;

Path_MainFolder = 'D:\Google\我的雲端硬碟\學術｜研究與論文\論文著作\MLE Pricing Kernel';
Path_Output      = fullfile(Path_MainFolder, 'Code', '11  Output');
Path_Output_Plot = fullfile(Path_MainFolder, 'Code', '11  Output - Plot');

%% Plot setting

set(groot, 'defaultAxesFontName', 'Times New Roman');
set(groot, 'defaultAxesFontSize', 12);
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultLineMarkerFaceColor','auto');


%% Plot M curve

% -------- User-specified (alpha, beta) by L --------
alpha_val = 1.00;
beta_val  = 0.90;

targets = [
    struct('L', 1, 'alpha', alpha_val, 'beta', beta_val);
    struct('L', 2, 'alpha', alpha_val, 'beta', beta_val);
    struct('L', 3, 'alpha', alpha_val, 'beta', beta_val)
];

select_rows = cell(3,1);
tol = 1e-9;
files = dir(fullfile(Path_Output, 'MLE_gamma_L_*.mat'));

for k = 1:numel(targets)
    L_i   = targets(k).L;
    a_tar = targets(k).alpha;
    b_tar = targets(k).beta;

    chosen = [];
    chosen_file = '';

    for f = 1:numel(files)
        Sfile = fullfile(Path_Output, files(f).name);
        S     = load(Sfile, 'L','alpha','beta','gamma_hat');

        if isfield(S,'L') && isfield(S,'alpha') && isfield(S,'beta')
            if S.L==L_i && abs(S.alpha-a_tar)<=tol && abs(S.beta-b_tar)<=tol
                chosen = S;
                chosen_file = files(f).name;
                break
            end
        end
    end
    if isempty(chosen)
        error('Exact file not found.');
    end
    select_rows{L_i} = table(L_i, a_tar, b_tar, string(chosen_file), ...
        'VariableNames', {'L','alpha','beta','file'});
end

% Define R_axis
R_axis = linspace(0.3, 3.0, 10000)'; 
logR   = log(R_axis);

% plot
figure('Position',[50 80 500 400]);
layout = tiledlayout(1, 1, 'TileSpacing', 'Compact', 'Padding', 'None');
nexttile;
hold on;

colors = get(groot, 'defaultAxesColorOrder');

for L_i = 3:3
    row = select_rows{L_i};

    Sfile = fullfile(Path_Output, row.file);
    data  = load(Sfile, 'gamma_hat', 'L', 'alpha', 'beta');
    gamma_hat = data.gamma_hat;
    L_current = data.L;
    fprintf('\nL = %d, alpha = %.2f, beta = %.2f, gamma_hat = %s\n', ...
        row.L, row.alpha, row.beta, mat2str(data.gamma_hat, 4));

    poly_sum = zeros(size(logR));
    for l = 1:L_current
        poly_sum = poly_sum + gamma_hat(l) * (logR.^l);
    end
    M_curve = exp(-poly_sum);

    plot(R_axis, M_curve, 'LineWidth', 1.8, ...
        'Color', colors(L_i,:), ...
        'DisplayName', sprintf('$L=%d$', row.L));
end

hold off;
xlabel('$R$');
ylabel('$M(R)$');
legend('show','Location','northeast','Box','off');
grid on;
xlim([0.8 1.2]);
set(gca,'LooseInset',get(gca,'TightInset'));

% Output
out_png = sprintf('plot_M_curve_alpha_%.2f_beta_%.2f.png', alpha_val, beta_val);
saveas(gcf, fullfile(Path_Output_Plot, out_png));


%% Compute risk preference indices and plot (given L, alpha, beta)

clc
colors = get(groot, 'defaultAxesColorOrder');

% -------- Parameter sets (six cases) --------
param_list = {
    struct('L',1,'alpha',0.95,'beta',0.90)
    struct('L',2,'alpha',0.95,'beta',0.90)
    struct('L',3,'alpha',1.00,'beta',0.90)
    struct('L',1,'alpha',1.00,'beta',1.00)
    struct('L',2,'alpha',1.00,'beta',1.00)
    struct('L',3,'alpha',1.00,'beta',1.00)
};

measures = {'ARA','RRA','AP','RP','AT','RT'};

x_min = 0.3;
x_max = 3.0;
num_points = 5000;

R_axis = linspace(x_min, x_max, num_points)';
x      = log(R_axis);

file_cache = containers.Map();
files = dir(fullfile(Path_Output, 'MLE_gamma_L_*.mat'));

% ============================================================
%  Main Loop
% ============================================================
% for idx = 1:length(param_list)
for idx = 3
    L_target     = param_list{idx}.L;
    alpha_target = param_list{idx}.alpha;
    beta_target  = param_list{idx}.beta;
    
    % ----- Find and Load file -----
    chosen_file = '';
    for f = 1:numel(files)
        Sfile = fullfile(Path_Output, files(f).name);
        S     = load(Sfile, 'L','alpha','beta','gamma_hat');
        if S.L == L_target && S.alpha == alpha_target && S.beta == beta_target
            chosen_file = files(f).name;
            gamma_hat = S.gamma_hat;
            break;
        end
    end
    if isempty(chosen_file)
        error('MLE file not found for L=%d alpha=%.2f beta=%.2f', ...
              L_target, alpha_target, beta_target);
    end
    fprintf('Processing: %s (L=%d)\n', chosen_file, L_target);
    
    % ----- Analytical Calculation -----
    % P(x) and its derivatives w.r.t x (where x = ln R)
    P              = zeros(size(x)); % P(x)
    P_prime        = zeros(size(x)); % P'(x)
    P_double_prime = zeros(size(x)); % P''(x)
    P_triple_prime = zeros(size(x)); % P'''(x)
    
    for l = 1:L_target
        % P(x) = sum( gamma_l * x^l )
        P = P + gamma_hat(l) * (x.^l);
        
        % P'(x)
        P_prime = P_prime + l * gamma_hat(l) * (x.^(l-1));
        
        if l >= 2
            % P''(x)
            P_double_prime = P_double_prime + l * (l-1) * gamma_hat(l) * (x.^(l-2));
        end
        
        if l >= 3
            % P'''(x)
            P_triple_prime = P_triple_prime + l * (l-1) * (l-2) * gamma_hat(l) * (x.^(l-3));
        end
    end
    
    % The SDF level: M(R) = exp(-P(ln R))
    % (ignoring constant exp(delta_t) as it cancels in ratios)
    M_val = exp(-P); 
    
    % --- Chain Rule Derivation for M derivatives ---
    % Let y = ln M = -P(x)
    % Derivatives of y w.r.t x:
    y_x   = -P_prime;
    y_xx  = -P_double_prime;
    y_xxx = -P_triple_prime;
    
    % Derivatives of y w.r.t R (using Chain Rule):
    % y' = dy/dR = (dy/dx) * (1/R)
    y1 = y_x ./ R_axis;
    
    % y'' = d^2y/dR^2
    y2 = (y_xx - y_x) ./ (R_axis.^2);
    
    % y''' = d^3y/dR^3
    y3 = (y_xxx - 3*y_xx + 2*y_x) ./ (R_axis.^3);
    
    % Calculate Absolute Derivatives M1, M2, M3
    % M1 = M'
    M1 = M_val .* y1;
    
    % M2 = M''
    M2 = M_val .* (y2 + y1.^2);
    
    % M3 = M'''
    M3 = M_val .* (y3 + 3*y1.*y2 + y1.^3);
    
    % --- Risk Indices (Computed directly from M, M1, M2, M3) ---
    % 1. ARA = -M' / M
    risk_pref.ARA = -M1 ./ M_val;
    
    % 2. RRA = R * ARA
    risk_pref.RRA = R_axis .* risk_pref.ARA;
    
    % 3. AP = -M'' / M'
    risk_pref.AP  = -M2 ./ M1;
    
    % 4. RP = R * AP
    risk_pref.RP  = R_axis .* risk_pref.AP;
    
    % 5. AT = -M''' / M''
    risk_pref.AT  = -M3 ./ M2;
    
    % 6. RT = R * AT
    risk_pref.RT  = R_axis .* risk_pref.AT;
    
    % ============================================================
    %  Output Data Table (Filtered & Formatted)
    % ============================================================
    T_full = table(R_axis, M_val, M1, M2, M3, ...
                   risk_pref.ARA, risk_pref.RRA, ...
                   risk_pref.AP, risk_pref.RP, ...
                   risk_pref.AT, risk_pref.RT, ...
                   'VariableNames', {'R', 'M', 'M1', 'M2', 'M3', ...
                                     'ARA', 'RRA', 'AP', 'RP', 'AT', 'RT'});
    
    mask = (R_axis >= 1.16) & (R_axis <= 1.19);
    T_out = T_full(mask, :);
    
    varNames = T_out.Properties.VariableNames;
    for v = 1:numel(varNames)
        T_out.(varNames{v}) = round(T_out.(varNames{v}), 8);
    end
    
    alpha_str = sprintf('%.2f', alpha_target);
    beta_str  = sprintf('%.2f', beta_target);
    csv_filename = sprintf('RiskTable_L_%d_alpha_%s_beta_%s.csv', L_target, alpha_str, beta_str);
    full_csv_path = fullfile(Path_Output, csv_filename);
    
    writetable(T_out, full_csv_path);
    fprintf('Saved filtered table: %s (Rows: %d)\n', csv_filename, height(T_out));
    
    % ============================================================
    %  Plotting (Using full data, NOT filtered data)
    % ============================================================
    for m = 1:length(measures)
        key    = measures{m};
        Y_plot = risk_pref.(key);
        
        fig = figure('Position',[50 80 450 400], 'Visible', 'off');
        layout = tiledlayout(1, 1, 'TileSpacing', 'Compact', 'Padding', 'None');
        nexttile;
        hold on;
        
        plot(R_axis, Y_plot, 'LineWidth', 1.5, ...
            'Color', colors(L_target,:));
        
        grid on;
        hold off;
        xlabel('$R$', 'Interpreter', 'latex');
        ylabel(key);
        
        axis tight; 
        switch key
            case 'ARA', y_limits = [0, 4.5]; 
            case 'RRA', y_limits = [0, 3.6];
            case 'AP',  y_limits = [-10, 100];
            case 'RP',  y_limits = [-10, 100];
            case 'AT',  y_limits = [-10, 100];
            case 'RT',  y_limits = [-10, 100];
            otherwise,  axis tight; y_limits = ylim;
        end
        ylim(y_limits);
        xlim([0.8 1.2]);
        
        set(gca,'LooseInset',get(gca,'TightInset'));
        
        out_png = sprintf('%s_L_%d_alpha_%s_beta_%s.png', ...
                          key, L_target, alpha_str, beta_str);
        saveas(fig, fullfile(Path_Output_Plot, out_png));
        close(fig);
    end
end
fprintf('\nAll plots and tables completed.\n');