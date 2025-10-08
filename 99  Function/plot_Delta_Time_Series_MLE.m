function [fig, date_dt] = plot_Delta_Time_Series_MLE(date_fields, delta_all, b_list, varargin)

    % ----- parse options -----
    ip = inputParser;
    addParameter(ip, 'YLim', [], @(x) isempty(x) || (isnumeric(x) && numel(x)==2));
    addParameter(ip, 'Path_Output', pwd, @ischar);
    addParameter(ip, 'FileName', 'Delta_t_Timeseries_MLE.png', @ischar);
    addParameter(ip, 'FigurePosition', [100 100 1500 500], @(x)isnumeric(x) && numel(x)==4);
    addParameter(ip, 'DateFormat', 'yyyyMMdd', @ischar);
    addParameter(ip, 'LineWidth', 1.2, @(x)isnumeric(x)&&isscalar(x));
    addParameter(ip, 'YLabel', '$\delta_t$', @ischar);
    addParameter(ip, 'YLabelRotation', 0, @(x)isnumeric(x)&&isscalar(x));
    parse(ip, varargin{:});
    opt = ip.Results;

    % ----- basic checks -----
    if numel(delta_all) ~= numel(b_list)
        error('delta_all (size %d) 必須與 b_list (size %d) 對應。', numel(delta_all), numel(b_list));
    end
    T = numel(date_fields);
    for j = 1:numel(b_list)
        if numel(delta_all{j}) ~= T
            error('delta_all{%d} 長度 (%d) 必須等於期數 T (%d)。', j, numel(delta_all{j}), T);
        end
    end

    % ----- time axis -----
    date_dt = datetime(date_fields, 'InputFormat', opt.DateFormat);

    % ----- plot -----
    fig = figure; set(fig, 'Position', opt.FigurePosition);
    tiledlayout(1, numel(b_list), 'TileSpacing','Compact','Padding','None');

    for j = 1:numel(b_list)
        nexttile;
        plot(date_dt, delta_all{j}, 'LineWidth', opt.LineWidth);
        grid on; box on;
        set(gca, 'FontName','Times New Roman');

        if ~isempty(opt.YLim), ylim(opt.YLim); end

        title(sprintf('$b$ = %d', b_list(j)), 'Interpreter','latex');
        if j == 1
            ylabel(opt.YLabel, 'Interpreter','latex', 'Rotation', opt.YLabelRotation);
        end
    end

    % ----- save -----
    if ~exist(opt.Path_Output, 'dir'), mkdir(opt.Path_Output); end
    out_png = fullfile(opt.Path_Output, opt.FileName);
    saveas(fig, out_png);
    disp(['Saved figure: ', out_png]);
end
