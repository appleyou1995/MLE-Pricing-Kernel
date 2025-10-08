function [fig, M_avg_all] = plot_Average_Pricing_Kernel( ...
    Smooth_AllR, M_all, b_list, R_base, varargin)

    % ---- parse options ----
    ip = inputParser;
    addParameter(ip, 'XLim', [0.8 1.2], @(x)isnumeric(x)&&numel(x)==2);
    addParameter(ip, 'YLim', [0.994 1.01], @(x)isnumeric(x)&&numel(x)==2);
    addParameter(ip, 'Path_Output', pwd, @ischar);
    addParameter(ip, 'FileName', 'Average_Pricing_Kernel_MLE.png', @ischar);
    addParameter(ip, 'FigurePosition', [120 120 680 480], @(x)isnumeric(x)&&numel(x)==4);
    addParameter(ip, 'LineWidth', 1.5, @(x)isnumeric(x)&&isscalar(x));
    addParameter(ip, 'LegendLocation', 'northeast', @ischar);
    addParameter(ip, 'Title', 'Average Pricing Kernel (MLE)', @ischar);
    addParameter(ip, 'XLabel', 'Gross Return $R$', @ischar);
    addParameter(ip, 'YLabel', '$\mathrm{E}(M)$', @ischar);
    addParameter(ip, 'InterpMethod', 'linear', @ischar);
    addParameter(ip, 'ExtrapVal', NaN, @(x)isnumeric(x)&&isscalar(x));
    parse(ip, varargin{:});
    opt = ip.Results;

    % ---- basics & checks ----
    date_fields = Smooth_AllR.Properties.VariableNames;
    T = numel(date_fields);
    J = numel(b_list);
    if numel(M_all) ~= J
        error('M_all length (%d) must match b_list length (%d).', numel(M_all), J);
    end
    R_base = R_base(:)';
    N = numel(R_base);

    % ---- compute Average M across all periods ----
    M_avg_all = NaN(J, N);
    for j = 1:J
        M_vec = M_all{j};
        if size(M_vec,1) ~= T
            error('M_all{%d} has %d periods, but Smooth_AllR has %d.', j, size(M_vec,1), T);
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
    hold on;
    for j = 1:J
        plot(R_base, M_avg_all(j,:), 'LineWidth', opt.LineWidth, ...
             'DisplayName', sprintf('$b = %d$', b_list(j)));
    end
    hold off;

    grid on; box on;
    xlim(opt.XLim); ylim(opt.YLim);
    xlabel(opt.XLabel, 'Interpreter','latex');
    ylabel(opt.YLabel, 'Interpreter','latex');
    title(opt.Title, 'Interpreter','latex');
    legend('Location', opt.LegendLocation, 'Box','off', 'Interpreter','latex');

    % ---- save ----
    if ~exist(opt.Path_Output, 'dir'), mkdir(opt.Path_Output); end
    saveas(fig, fullfile(opt.Path_Output, opt.FileName));
end
