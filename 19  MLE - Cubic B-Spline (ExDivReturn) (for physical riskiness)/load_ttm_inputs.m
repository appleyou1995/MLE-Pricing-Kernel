function Data = load_ttm_inputs(Path_Data, Path_RND, Risk_Free_Rate_All, ...
    Target_TTM, Sample_Start, Sample_End, Expected_Grid_Size)
%LOAD_TTM_INPUTS Load one TTM using its own Hsieh contract-date sample.

    target_file = fullfile(Path_Data, sprintf('TTM_%d.csv', Target_TTM));
    Target_Table = readtable(target_file);
    require_variables(Target_Table, {'date', 'exdate'}, target_file);

    target_dates = double(Target_Table.date);
    target_exdates = double(Target_Table.exdate);
    in_sample = target_dates >= Sample_Start & target_dates <= Sample_End;
    target_dates = target_dates(in_sample);
    target_exdates = target_exdates(in_sample);

    [target_dates, target_sort_idx] = sort(target_dates);
    target_exdates = target_exdates(target_sort_idx);

    if isempty(target_dates)
        error('TTM = %d has no target dates in the requested sample.', Target_TTM);
    end
    if numel(unique(target_dates)) ~= numel(target_dates)
        error('Duplicate target dates found in %s.', target_file);
    end
    target_months = floor(target_dates ./ 100);
    if numel(unique(target_months)) ~= numel(target_months)
        error('TTM = %d has more than one target date in a calendar month.', ...
            Target_TTM);
    end

    required_years = unique(floor(target_dates ./ 10000));
    Smooth_AllR = table();
    Smooth_AllR_RND = table();

    for calendar_year = required_years(:)'
        input_file = fullfile(Path_RND, sprintf( ...
            'TTM_%d_RND_Tables_%d.mat', Target_TTM, calendar_year));
        if ~isfile(input_file)
            error('Required RND file is missing: %s', input_file);
        end

        S = load(input_file, 'Table_Smooth_AllR', 'Table_Smooth_AllR_RND');
        if ~isfield(S, 'Table_Smooth_AllR') || ...
                ~isfield(S, 'Table_Smooth_AllR_RND')
            error('RND file lacks required tables: %s', input_file);
        end
        if width(S.Table_Smooth_AllR) ~= width(S.Table_Smooth_AllR_RND)
            error('R and RND column counts differ in %s.', input_file);
        end

        Smooth_AllR = [Smooth_AllR, S.Table_Smooth_AllR];             %#ok<AGROW>
        Smooth_AllR_RND = [Smooth_AllR_RND, S.Table_Smooth_AllR_RND]; %#ok<AGROW>
    end

    rnd_var_names = string(Smooth_AllR.Properties.VariableNames);
    rnd_dates = NaN(numel(rnd_var_names), 1);
    for j = 1:numel(rnd_var_names)
        date_token = regexp(rnd_var_names(j), '\d{8}', 'match', 'once');
        if ~isempty(date_token)
            rnd_dates(j) = str2double(date_token);
        end
    end
    if any(~isfinite(rnd_dates))
        bad_names = strjoin(rnd_var_names(~isfinite(rnd_dates)), ', ');
        error('Cannot parse quote dates from RND variable names: %s', bad_names);
    end
    if numel(unique(rnd_dates)) ~= numel(rnd_dates)
        error('Duplicate RND quote dates found for TTM = %d.', Target_TTM);
    end

    in_sample = rnd_dates >= Sample_Start & rnd_dates <= Sample_End;
    Smooth_AllR = Smooth_AllR(:, in_sample);
    Smooth_AllR_RND = Smooth_AllR_RND(:, in_sample);
    rnd_dates = rnd_dates(in_sample);

    [rnd_dates, sort_idx] = sort(rnd_dates);
    Smooth_AllR = Smooth_AllR(:, sort_idx);
    Smooth_AllR_RND = Smooth_AllR_RND(:, sort_idx);

    missing_rnd_dates = setdiff(target_dates, rnd_dates);
    extra_rnd_dates = setdiff(rnd_dates, target_dates);
    if ~isempty(missing_rnd_dates) || ~isempty(extra_rnd_dates)
        error(['TTM = %d RND dates do not exactly match TTM_%d.csv. ' ...
            'Missing RND dates: %s; unexpected RND dates: %s'], ...
            Target_TTM, Target_TTM, preview_dates(missing_rnd_dates), ...
            preview_dates(extra_rnd_dates));
    end
    if ~isequal(rnd_dates(:), target_dates(:))
        error('TTM = %d date ordering could not be aligned.', Target_TTM);
    end

    realized_file = fullfile(Path_Data, sprintf( ...
        'Realized_Return_TTM_%d.csv', Target_TTM));
    Realized_All = readtable(realized_file);
    require_variables(Realized_All, {'date', 'realized_ret'}, realized_file);
    realized_dates = double(Realized_All.date);
    if numel(unique(realized_dates)) ~= numel(realized_dates)
        error('Duplicate dates found in %s.', realized_file);
    end

    rf_date_col = sprintf('date_%d', Target_TTM);
    rf_rate_col = sprintf('rf_gross_TTM%d', Target_TTM);
    require_variables(Risk_Free_Rate_All, {rf_date_col, rf_rate_col}, ...
        'Risk_Free_GrossFactor_ByTargetTTM.csv');
    rf_dates_all = double(Risk_Free_Rate_All.(rf_date_col));
    rf_values_all = double(Risk_Free_Rate_All.(rf_rate_col));
    rf_valid = isfinite(rf_dates_all) & isfinite(rf_values_all);
    rf_dates_all = rf_dates_all(rf_valid);
    rf_values_all = rf_values_all(rf_valid);
    if numel(unique(rf_dates_all)) ~= numel(rf_dates_all)
        error('Duplicate risk-free dates found for TTM = %d.', Target_TTM);
    end

    assert_all_dates_present(rnd_dates, realized_dates, ...
        'realized return', Target_TTM);
    assert_all_dates_present(rnd_dates, rf_dates_all, ...
        'risk-free rate', Target_TTM);

    [~, idx_realized] = ismember(rnd_dates, realized_dates);
    [~, idx_rf] = ismember(rnd_dates, rf_dates_all);

    Realized_Return = Realized_All(idx_realized, :);
    Risk_Free_Rate = rf_values_all(idx_rf);
    Expiration_Dates = target_exdates;

    quote_dt = datetime(string(rnd_dates), 'InputFormat', 'yyyyMMdd');
    expiration_dt = datetime(string(Expiration_Dates), ...
        'InputFormat', 'yyyyMMdd');
    Actual_TTM_Days = days(expiration_dt - quote_dt);

    Grid_Lengths = zeros(numel(rnd_dates), 1);
    Raw_Min_R = Inf;
    Raw_Max_R = -Inf;
    fields = Smooth_AllR.Properties.VariableNames;
    rnd_fields = Smooth_AllR_RND.Properties.VariableNames;

    for j = 1:numel(fields)
        R_axis = double(Smooth_AllR.(fields{j})(:));
        rnd_pdf = double(Smooth_AllR_RND.(rnd_fields{j})(:));
        Grid_Lengths(j) = numel(R_axis);

        if numel(rnd_pdf) ~= numel(R_axis)
            error('R and RND lengths differ on quote date %08d.', rnd_dates(j));
        end
        if Expected_Grid_Size > 0 && numel(R_axis) ~= Expected_Grid_Size
            error('Quote date %08d has %d grid points; expected %d.', ...
                rnd_dates(j), numel(R_axis), Expected_Grid_Size);
        end
        if any(~isfinite(R_axis)) || any(~isfinite(rnd_pdf))
            error('Non-finite R or RND value found on quote date %08d.', ...
                rnd_dates(j));
        end
        if any(diff(R_axis) <= 0)
            error('R grid is not strictly increasing on quote date %08d.', ...
                rnd_dates(j));
        end
        if any(rnd_pdf < 0) || trapz(R_axis, rnd_pdf) <= 0
            error('Invalid RND density found on quote date %08d.', rnd_dates(j));
        end

        Raw_Min_R = min(Raw_Min_R, R_axis(1));
        Raw_Max_R = max(Raw_Max_R, R_axis(end));
    end

    range_R = Raw_Max_R - Raw_Min_R;
    Global_Min_R = max(eps, Raw_Min_R - 0.10 * range_R);
    Global_Max_R = Raw_Max_R + 0.10 * range_R;

    Data = struct();
    Data.Smooth_AllR = Smooth_AllR;
    Data.Smooth_AllR_RND = Smooth_AllR_RND;
    Data.Realized_Return = Realized_Return;
    Data.Risk_Free_Rate = Risk_Free_Rate(:);
    Data.Quote_Dates = rnd_dates(:);
    Data.Expiration_Dates = Expiration_Dates(:);
    Data.Actual_TTM_Days = Actual_TTM_Days(:);
    Data.Grid_Lengths = Grid_Lengths;
    Data.Raw_Min_R = Raw_Min_R;
    Data.Raw_Max_R = Raw_Max_R;
    Data.Global_Min_R = Global_Min_R;
    Data.Global_Max_R = Global_Max_R;
end


function require_variables(T, required_names, source_name)
    available = T.Properties.VariableNames;
    missing = required_names(~ismember(required_names, available));
    if ~isempty(missing)
        error('%s is missing variables: %s', source_name, ...
            strjoin(missing, ', '));
    end
end


function assert_all_dates_present(master_dates, available_dates, label, Target_TTM)
    missing = setdiff(master_dates, available_dates);
    if ~isempty(missing)
        error('TTM = %d is missing %s data for quote date(s): %s', ...
            Target_TTM, label, preview_dates(missing));
    end
end


function output = preview_dates(date_values)
    if isempty(date_values)
        output = 'none';
        return
    end
    date_values = date_values(:);
    shown = date_values(1:min(6, numel(date_values)));
    output = strtrim(sprintf('%08d ', shown));
    if numel(date_values) > numel(shown)
        output = [output ' ...'];
    end
end
