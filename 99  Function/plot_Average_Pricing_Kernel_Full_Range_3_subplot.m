function [fig, M_avg_all] = plot_Average_Pricing_Kernel_Full_Range_3_subplot( ...
    Smooth_AllR, M_all, param_count_list, R_base, varargin)

    % ---- parse options ----
    ip = inputParser;
    addParameter(ip, 'XLim', [0 3], @(x)isnumeric(x)&&numel(x)==2);
    addParameter(ip, 'YLim', [], @(x) isempty(x) || (isnumeric(x)&&numel(x)==2));
    addParameter(ip, 'InterpMethod', 'linear', @ischar);
    addParameter(ip, 'ExtrapVal', NaN, @(x)isnumeric(x)&&isscalar(x));
    addParameter(ip, 'Path_Output', pwd, @ischar);
    addParameter(ip, 'FileName', 'Average_Pricing_Kernel_Full_Range.png', @ischar);
    addParameter(ip, 'FigurePosition', [100 100 1200 400], @(x)isnumeric(x)&&numel(x)==4);
    addParameter(ip, 'TileSpacing', 'Compact', @ischar);
    addParameter(ip, 'Padding', 'None', @ischar);
    addParameter(ip, 'LineWidth', 1.5, @(x)isnumeric(x)&&isscalar(x));
    addParameter(ip, 'FontName', 'Times New Roman', @ischar);
    addParameter(ip, 'Interpreter', 'latex', @ischar);
    parse(ip, varargin{:});
    opt = ip.Results;

    % ---- basics & checks ----
    date_fields = Smooth_AllR.Properties.VariableNames;
    T = numel(date_fields);
    J = numel(param_count_list);
    if numel(M_all) ~= J
        error('M_all 的長度 (%d) 必須等於 param_count_list 的長度 (%d)。', numel(M_all), J);
    end
    R_base = R_base(:)';
    N = numel(R_base);

    % ---- compute average M for each b ----
    M_avg_all = NaN(J, N);

    for j = 1:J
        M_vec = M_all{j};
        if size(M_vec,1) ~= T
            error('M_all{%d} 的期數 (%d) 與 Smooth_AllR 期數 (%d) 不符。', j, size(M_vec,1), T);
        end

        M_interp = NaN(T, N);
        for t = 1:T
            R_t = Smooth_AllR.(date_fields{t});
            M_t = M_vec(t,:);
            M_interp(t,:) = interp1(R_t, M_t, R_base, opt.InterpMethod, opt.ExtrapVal);
        end
        M_avg_all(j,:) = mean(M_interp, 1, 'omitnan');
    end

    % ---- plot ----
    fig = figure; set(fig, 'Position', opt.FigurePosition);
    tiledlayout(1, J, 'TileSpacing', opt.TileSpacing, 'Padding', opt.Padding);

    for j = 1:J
        nexttile;
        plot(R_base, M_avg_all(j,:), 'LineWidth', opt.LineWidth);
        set(gca, 'FontName', opt.FontName);
        grid on; xlim(opt.XLim);
        if ~isempty(opt.YLim), ylim(opt.YLim); end

        xlabel('Gross Return', 'FontName', opt.FontName);
        if j == 1
            ylabel('$\mathrm{E}(M)$', 'Interpreter', opt.Interpreter);
        end
        title(sprintf('$L=%d$', param_count_list(j)), 'FontName', opt.FontName, 'Interpreter', opt.Interpreter);
    end

    % ---- save ----
    if ~exist(opt.Path_Output, 'dir'), mkdir(opt.Path_Output); end
    saveas(fig, fullfile(opt.Path_Output, opt.FileName));
end
