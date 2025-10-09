function [Tcmd, LL, GammaMat, max_L_list, gammas] = build_MLE_summary_Boswijk(Path_Output)

    files = dir(fullfile(Path_Output, 'MLE_gamma_max_L_*.mat'));
    if isempty(files)
        error('%s can not find MLE_gamma_max_L_*.mat。', Path_Output);
    end
    
    max_L_vals = nan(1, numel(files));
    keep = false(1, numel(files));
    
    for k = 1:numel(files)
        fname = files(k).name;
    
        L_value = regexp(fname, '(?<=_L_)(\d+)', 'match', 'once');
    
        if ~isempty(L_value)
            max_L_vals(k) = str2double(L_value);
            keep(k) = true;
        else
            warning('Skip file (cannot parse L index): %s', fname);
        end
    end
    
    files = files(keep);
    max_L_vals = max_L_vals(keep);
    
    if isempty(files)
        error('No valid files with an L index found in: %s', Path_Output);
    end
    
    [max_L_list, order] = sort(max_L_vals);
    files = files(order);

    J = numel(files);
    LL = nan(1, J);
    gammas = cell(1, J);
    max_L = max(max_L_list);

    for j = 1:J
        S = load(fullfile(Path_Output, files(j).name));
        if ~isfield(S, 'log_lik') || ~isfield(S, 'gamma_hat')
            error('missing log_lik or gamma_hat：%s', files(j).name);
        end
        LL(j) = S.log_lik;
        gammas{j} = S.gamma_hat(:);
    end

    GammaMat = zeros(max_L, J);
    for j = 1:J
        Lj = max_L_list(j);
        g = gammas{j}(:);
        g = g(1:Lj);
        GammaMat(1:Lj, j) = g;
    end

    fmt = @(A) A .* (abs(A) >= 5e-5);
    LL_disp = fmt(LL);
    Gamma_disp = fmt(GammaMat);

    rowNames = ["Log-Likelihood"; compose('\\gamma_{%d}', (1:max_L)')];
    varNames = compose('max_L_%d', max_L_list);
    Tcmd = array2table([LL_disp; Gamma_disp], 'RowNames', rowNames, 'VariableNames', varNames);

    disp('===== MLE summary =====');
    disp(Tcmd);
end
