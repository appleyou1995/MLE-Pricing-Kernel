function [fig, M_avg_all, band_lo_all, band_hi_all] = ...
    plot_Average_Pricing_Kernel_Extreme_Interval_3_subplot( ...
        Smooth_AllR, M_all, b_list, R_base, varargin)

    % ---- options ----
    ip = inputParser;
    addParameter(ip, 'XLim', [0.8 1.2], @(x)isnumeric(x)&&numel(x)==2);
    addParameter(ip, 'YLim', [], @(x) isempty(x) || (isnumeric(x)&&numel(x)==2));
    addParameter(ip, 'InterpMethod', 'linear', @ischar);
    addParameter(ip, 'ExtrapVal', NaN, @(x)isnumeric(x)&&isscalar(x));
    addParameter(ip, 'Path_Output', pwd, @ischar);
    addParameter(ip, 'FileName', 'Average_Pricing_Kernel_Interval.png', @ischar);
    addParameter(ip, 'FigurePosition', [100 100 1200 400], @(x)isnumeric(x)&&numel(x)==4);
    addParameter(ip, 'TileSpacing', 'Compact', @ischar);
    addParameter(ip, 'Padding', 'None', @ischar);
    addParameter(ip, 'LineWidth', 1.5, @(x)isnumeric(x)&&isscalar(x));
    addParameter(ip, 'FontName', 'Times New Roman', @ischar);
    addParameter(ip, 'Interpreter', 'latex', @ischar);
    addParameter(ip, 'LowerPct', 5,  @(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=100);
    addParameter(ip, 'UpperPct', 95, @(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=100);
    addParameter(ip, 'BandAlpha', 0.1, @(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=1);
    parse(ip, varargin{:});
    opt = ip.Results;

    % ---- checks ----
    date_fields = Smooth_AllR.Properties.VariableNames;
    T = numel(date_fields);
    J = numel(b_list);
    if numel(M_all) ~= J
        error('M_all 的長度 (%d) 必須等於 b_list 的長度 (%d)。', numel(M_all), J);
    end
    R_base = R_base(:)';  N = numel(R_base);

    % ---- compute: mean & percentile bands across periods ----
    M_avg_all   = NaN(J, N);
    band_lo_all = NaN(J, N);
    band_hi_all = NaN(J, N);

    for j = 1:J
        M_vec = M_all{j};                 % T×N_t
        if size(M_vec,1) ~= T
            error('M_all{%d} 的期數 (%d) 與 Smooth_AllR 期數 (%d) 不符。', j, size(M_vec,1), T);
        end

        M_interp = NaN(T, N);
        for t = 1:T
            R_t = Smooth_AllR.(date_fields{t});
            M_t = M_vec(t,:);
            M_interp(t,:) = interp1(R_t, M_t, R_base, opt.InterpMethod, opt.ExtrapVal);
        end

        % 跨期平均
        M_avg_all(j,:) = mean(M_interp, 1, 'omitnan');

        % 百分位（逐列：沿時間维度取百分位），忽略 NaN
        [band_lo_all(j,:), band_hi_all(j,:)] = prctile_cols_omitnan( ...
            M_interp, opt.LowerPct, opt.UpperPct);
    end

    % ---- plot ----
    fig = figure; set(fig, 'Position', opt.FigurePosition);
    tiledlayout(1, J, 'TileSpacing', opt.TileSpacing, 'Padding', opt.Padding);

    for j = 1:J
        nexttile; hold on;

        % 先畫帶狀區（用同色半透明）
        lo = band_lo_all(j,:);  hi = band_hi_all(j,:);
        mask = isfinite(lo) & isfinite(hi) & (hi >= lo);
        if any(mask)
            % 單段填色（若 mask 非連續，視覺上仍足夠）
            xpoly = [R_base(mask), fliplr(R_base(mask))];
            ypoly = [lo(mask),     fliplr(hi(mask))];
            p = patch(xpoly, ypoly, [0 0 0], 'EdgeColor','none', 'FaceAlpha', opt.BandAlpha);
        end

        % 畫平均線
        hMean = plot(R_base, M_avg_all(j,:), 'LineWidth', opt.LineWidth, ...
                     'DisplayName', sprintf('$b = %d$', b_list(j)));

        % 以平均線的顏色畫上下百分位虛線
        col = get(hMean, 'Color');
        if exist('p','var') && isgraphics(p), set(p, 'FaceColor', col); end
        plot(R_base, lo, '--', 'Color', col, 'LineWidth', 1.0, 'HandleVisibility','off');
        plot(R_base, hi, '--', 'Color', col, 'LineWidth', 1.0, 'HandleVisibility','off');

        % 軸&標題
        set(gca, 'FontName', opt.FontName);
        grid on; xlim(opt.XLim);
        if ~isempty(opt.YLim), ylim(opt.YLim); end
        xlabel('Gross Return', 'FontName', opt.FontName);
        if j == 1
            ylabel('$\mathrm{E}(M)$', 'Interpreter', opt.Interpreter);
        end
        title(sprintf('$b=%d$', b_list(j)), 'FontName', opt.FontName, 'Interpreter', opt.Interpreter);

        hold off;
        clear p
    end

    % ---- save ----
    if ~exist(opt.Path_Output, 'dir'), mkdir(opt.Path_Output); end
    saveas(fig, fullfile(opt.Path_Output, opt.FileName));
end

% ===== helper: percentile-by-column with NaN omission =====
function [lo, hi] = prctile_cols_omitnan(X, pLo, pHi)
    N = size(X,2);
    lo = NaN(1,N);  hi = NaN(1,N);
    for k = 1:N
        col = X(:,k);
        col = col(isfinite(col));
        if ~isempty(col)
            lo(k) = prctile(col, pLo);
            hi(k) = prctile(col, pHi);
        end
    end
end
