function [InvE_M_all, bad_count, date_dt, fig] = plot_InvE_M_vs_Rf( ...
    Smooth_AllR, Smooth_AllR_RND, Risk_Free_Rate, M_all, b_list, varargin)

    % ---- parse options ----
    ip = inputParser;
    addParameter(ip, 'YLim', [], @(x) isempty(x) || (isnumeric(x)&&numel(x)==2));
    addParameter(ip, 'Path_Output', pwd, @ischar);
    addParameter(ip, 'FileName', 'InvE_M_vs_Rf_Timeseries_MLE.png', @ischar);
    addParameter(ip, 'FigurePosition', [100 100 1500 500], @(x)isnumeric(x)&&numel(x)==4);
    addParameter(ip, 'DateFormat', 'yyyyMMdd', @ischar);
    addParameter(ip, 'MaskFrac', 0.6, @(x)isnumeric(x)&&isscalar(x)&&x>0&&x<=1);
    parse(ip, varargin{:});
    opt = ip.Results;

    % ---- basic checks ----
    date_fields = Smooth_AllR.Properties.VariableNames;
    T = numel(date_fields);
    if numel(Risk_Free_Rate) ~= T
        error('Risk_Free_Rate length (%d) != num periods in Smooth_AllR (%d).', numel(Risk_Free_Rate), T);
    end
    if numel(M_all) ~= numel(b_list)
        error('M_all size (%d) must match b_list size (%d).', numel(M_all), numel(b_list));
    end
    for j = 1:numel(b_list)
        if size(M_all{j},1) ~= T
            error('M_all{%d} first dimension (%d) != T (%d).', j, size(M_all{j},1), T);
        end
    end

    % ---- compute 1 / E[M] per period ----
    date_dt = datetime(date_fields, 'InputFormat', opt.DateFormat);
    InvE_M_all = nan(T, numel(b_list));
    bad_count  = zeros(1, numel(b_list));

    for j = 1:numel(b_list)
        M_mat = M_all{j};   % T×N
        for t = 1:T
            R  = Smooth_AllR.(date_fields{t})(:);       % N×1
            fs = Smooth_AllR_RND.(date_fields{t})(:);   % N×1
            M  = M_mat(t, :)';                          % N×1
            Rf = Risk_Free_Rate(t);

            % 1) normalize f*_t
            Z = trapz(R, fs);
            if ~(Z > 0) || ~isfinite(Z), bad_count(j)=bad_count(j)+1; continue; end
            fs = fs ./ Z;

            % 2) valid mask: denom > 0 & finite
            denom = Rf .* M;
            mask  = isfinite(R) & isfinite(fs) & isfinite(M) & (denom > 0);
            if nnz(mask) < opt.MaskFrac * numel(R)
                bad_count(j) = bad_count(j) + 1; continue
            end

            % 3) compute f_t and E[M]
            f_phys = fs(mask) ./ denom(mask);
            EM_t   = trapz(R(mask), M(mask) .* f_phys);   % ≈ 1 / Rf
            if isfinite(EM_t) && EM_t > 0
                InvE_M_all(t, j) = 1.0 / EM_t;
            else
                bad_count(j) = bad_count(j) + 1;
            end
        end
    end

    fprintf('Invalid periods: ');
    fprintf('b=%d: %d  ', [b_list; bad_count]);
    fprintf('(of %d)\n', T);

    % ---- plot ----
    fig = figure; set(fig, 'Position', opt.FigurePosition);
    tiledlayout(1, numel(b_list), 'TileSpacing','Compact','Padding','None');

    for j = 1:numel(b_list)
        nexttile; hold on;
        plot(date_dt, InvE_M_all(:,j), 'LineWidth', 2, 'DisplayName','$1/E[M]$');
        plot(date_dt, Risk_Free_Rate,  '--', 'LineWidth', 2, 'DisplayName','$R_t^{f}$');
        hold off; grid on; box on;
        set(gca, 'FontName','Times New Roman');

        if ~isempty(opt.YLim)
            ylim(opt.YLim);
        else
            ys = [InvE_M_all(:,j); Risk_Free_Rate]; ys = ys(isfinite(ys));
            if ~isempty(ys)
                pad = 0.02 * (max(ys)-min(ys) + eps);
                ylim([min(ys)-pad, max(ys)+pad]);
            end
        end

        title(sprintf('$b$ = %d', b_list(j)), 'Interpreter','latex');
        if j == numel(b_list)
            legend('Location','best','Box','off','Interpreter','latex');
        end
    end

    % ---- save ----
    if ~exist(opt.Path_Output, 'dir'), mkdir(opt.Path_Output); end
    out_png = fullfile(opt.Path_Output, opt.FileName);
    saveas(fig, out_png);
    disp(['Saved figure: ', out_png]);
end
