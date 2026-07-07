function RF_Target = Function__RF_TTM(Data, Date_Target, TTM_Target)
    
    % Data:
    % [1. Date (YYYYMMDD) | 2. TTM (Days) | 3. RF (Annualized)]
    %
    % Date_Target:
    % YYYYMMDD numeric vector
    %
    % TTM_Target:
    % target TTM in days
    %
    % Output:
    % RF_Target:
    % interpolated annualized risk-free rate

    % Unique date-TTM pairs
    [AllPair, ~, Index_Pair] = unique([Date_Target(:), TTM_Target(:)], 'rows');

    RF_Pair = nan(size(AllPair, 1), 1);

    % Convert RF dates once
    Date_Data = datetime(string(Data(:, 1)), 'InputFormat', 'yyyyMMdd');

    for d = 1:size(AllPair, 1)

        targetDate = datetime(string(AllPair(d, 1)), 'InputFormat', 'yyyyMMdd');
        targetTTM  = AllPair(d, 2);

        found = false;

        % Search target date and up to 10 calendar days before
        for i = 0:10

            searchDate = targetDate - days(i);
            sameDate = Date_Data == searchDate;

            if any(sameDate)

                Data_TTM = Data(sameDate, 2);
                Data_RF  = Data(sameDate, 3);

                % Sort by TTM
                [Data_TTM, idx_sort] = sort(Data_TTM);
                Data_RF = Data_RF(idx_sort);

                % Average duplicated TTM if any
                [Data_TTM_U, ~, idx_group] = unique(Data_TTM);
                Data_RF_U = accumarray(idx_group, Data_RF, [], @mean);

                found = true;
                break

            end
        end

        if ~found
            warning('No RF data found for target date %s. Searched up to 10 days before.', ...
                    string(targetDate, 'yyyy-MM-dd'));
            continue
        end

        % Linear interpolation or extrapolation
        if targetTTM >= min(Data_TTM_U) && targetTTM <= max(Data_TTM_U)

            RF_Pair(d) = interp1(Data_TTM_U, Data_RF_U, targetTTM, 'linear');

        else

            RF_Pair(d) = interp1(Data_TTM_U, Data_RF_U, targetTTM, 'linear', 'extrap');

            warning('RF extrapolated for date %s, TTM = %.0f.', ...
                    string(targetDate, 'yyyy-MM-dd'), targetTTM);

        end

    end

    RF_Target = RF_Pair(Index_Pair);

end