function [Tcmd, LL, ThetaMat, b_list, thetas] = build_MLE_summary(Path_Output)

    files = dir(fullfile(Path_Output, 'MLE_theta_b*.mat'));
    if isempty(files)
        error('%s can not find MLE_theta_b*.mat。', Path_Output);
    end

    b_vals = nan(1, numel(files));
    for k = 1:numel(files)
        tok = regexp(files(k).name, '(?<=_b)(\d+)', 'match', 'once');
        if isempty(tok)
            error('The file name cannot be parsed b：%s', files(k).name);
        end
        b_vals(k) = str2double(tok);
    end
    [b_list, order] = sort(b_vals);
    files = files(order);

    J = numel(files);
    LL = nan(1, J);
    thetas = cell(1, J);
    max_b = max(b_list);

    for j = 1:J
        S = load(fullfile(Path_Output, files(j).name));
        if ~isfield(S, 'log_lik') || ~isfield(S, 'theta_hat')
            error('missing log_lik or theta_hat：%s', files(j).name);
        end
        LL(j) = S.log_lik;
        thetas{j} = S.theta_hat(:);
    end

    ThetaMat = zeros(max_b + 1, J);
    for j = 1:J
        bj = b_list(j);
        th = thetas{j};
        ThetaMat(1:bj+1, j) = th(1:bj+1);
        if bj < max_b
            ThetaMat(bj+2:end, j) = 0;
        end
    end

    fmt = @(A) A .* (abs(A) >= 5e-5);
    LL_disp = fmt(LL);
    Theta_disp = fmt(ThetaMat);

    rowNames = ["Log-Likelihood"; compose('\\theta_{%d}', (0:max_b)')];
    varNames = compose('b%d', b_list);
    Tcmd = array2table([LL_disp; Theta_disp], 'RowNames', rowNames, 'VariableNames', varNames);

    disp('===== MLE summary =====');
    disp(Tcmd);
end
