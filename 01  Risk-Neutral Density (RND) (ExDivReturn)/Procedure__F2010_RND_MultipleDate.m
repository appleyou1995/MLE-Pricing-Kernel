clear; clc

Path_MainFolder = 'D:\Google\我的雲端硬碟\學術｜研究與論文\論文著作\MLE Pricing Kernel';

Path_Data     = fullfile(Path_MainFolder, 'Code', '00  Output');
Path_Output   = fullfile(Path_MainFolder, 'Code', '01  Output');
Path_Data_inc = fullfile(Path_MainFolder, 'Data', 'IndexOptions1996202508_SP500', 'IV-Based');

% Specific Time-to-Maturity 
% Target_AllTTM = 30;
Target_AllTTM = 60;
% Target_AllTTM = 90;
% Target_AllTTM = 180;

Num_Grid = 30000;


%% Load Data

% Risk-Free Rate (Annualized)
% [1. Date (YYYYMMDD) | 2. TTM (Days) | 3. Risk-Free Rate (Annualized)]
FileName = 'RiskFreeRate19962025.txt';
Data_RF = load(fullfile(Path_Data, FileName));
clear FileName

% S&P 500 Dividend Yield
% [1. SecID | 2. Date (YYYYMMDD) | 3. Dividend Yield (Annualized)]
FileName = 'IndexDivYield19962025.txt';
Data_DY = load(fullfile(Path_Data, FileName));
clear FileName

% Remove SecID column
Data_DY(:, 1) = [];

% Convert dividend date to datenum for fallback search
Data_DY_DateNum = datenum(num2str(Data_DY(:, 1)), 'yyyymmdd');


%% Load Target Date

FileName = ['TTM_', num2str(Target_AllTTM), '.csv'];
Target_AllDate = readtable(fullfile(Path_Data, FileName));
clear FileName

% YYYYMMDD numeric
Target_AllDate.date   = str2double(string(Target_AllDate.date));
Target_AllDate.exdate = str2double(string(Target_AllDate.exdate));

% Compute calendar TTM
Target_Date_dt   = datetime(string(Target_AllDate.date),   'InputFormat', 'yyyyMMdd');
Target_ExDate_dt = datetime(string(Target_AllDate.exdate), 'InputFormat', 'yyyyMMdd');

Target_AllDate.TTM = days(Target_ExDate_dt - Target_Date_dt);

clear Target_Date_dt Target_ExDate_dt


%%  Test Mode: Select limited target dates

% Options:
% "single_year_1996" : 只跑 1996 年
% "special_dates"    : 只跑 Good Friday / Juneteenth / 近期測試日期
% "all_dates"        : 跑完整樣本

% RUN_MODE = "special_dates";
% RUN_MODE = "single_year_1996";
RUN_MODE = "all_dates";

switch RUN_MODE

    case "single_year_1996"

        Target_AllDate = Target_AllDate(Target_AllDate.date >= 19960101 & ...
                                        Target_AllDate.date <= 19961231, :);

        disp('RUN_MODE = single_year_1996');

        Target_AllDate_Display = table( ...
            string(compose('%08.0f', Target_AllDate.date)), ...
            string(compose('%08.0f', Target_AllDate.exdate)), ...
            Target_AllDate.TTM, ...
            'VariableNames', {'date', 'exdate', 'TTM'} ...
        );
        
        disp(Target_AllDate_Display);
        
        clear Target_AllDate_Display

    case "special_dates"

        % TTM_30 target dates:
        %
        % Normal first sample:
        % 19960117 -> 19960216
        %
        % Good Friday cases:
        % 20000321 -> 20000420
        % 20030318 -> 20030417
        % 20080219 -> 20080320
        % 20140318 -> 20140417
        % 20190319 -> 20190418
        % 20220315 -> 20220414
        % 20250318 -> 20250417
        %
        % Juneteenth quote-date shift:
        % 20240618 -> 20240719
        %
        % Recent normal case:
        % 20250716 -> 20250815

        Test_Date_List = [
            19960117
            20000321
            20030318
            20080219
            20140318
            20190319
            20220315
            20240618
            20250318
            20250716
        ];

        idx_test = ismember(Target_AllDate.date, Test_Date_List);

        Missing_Test_Dates = setdiff(Test_Date_List, Target_AllDate.date);

        if ~isempty(Missing_Test_Dates)
            warning('Some test dates are not in Target_AllDate: %s', ...
                    strjoin(string(Missing_Test_Dates'), ', '));
        end

        Target_AllDate = Target_AllDate(idx_test, :);

        disp('RUN_MODE = special_dates');

        Target_AllDate_Display = table( ...
            string(compose('%08.0f', Target_AllDate.date)), ...
            string(compose('%08.0f', Target_AllDate.exdate)), ...
            Target_AllDate.TTM, ...
            'VariableNames', {'date', 'exdate', 'TTM'} ...
        );
        
        disp(Target_AllDate_Display);
        
        clear Target_AllDate_Display

        clear Test_Date_List idx_test Missing_Test_Dates


    case "all_dates"

        disp('RUN_MODE = all_dates');

    otherwise

        error('Unknown RUN_MODE.');

end


%%  Target dates after test-mode filtering

Target_AllDate = sortrows(Target_AllDate, 'date');

Target_AllDate_date   = Target_AllDate.date;
Target_AllDate_exdate = Target_AllDate.exdate;

years = unique(floor(Target_AllDate_date / 10000));


%% Construct Risk-Neutral Density

Path_Data_01 = fullfile(Path_MainFolder, 'Code', '01  Risk-Neutral Density (RND)');
addpath(Path_Data_01);

Table_Diagnostics_All = table();

for y = 1:length(years)

    year_now = years(y);
    disp(['Processing year: ', num2str(year_now)]);
    
    Table_Smooth_AllK     = table();
    Table_Smooth_AllR     = table();
    Table_Smooth_AllR_RND = table();
    Table_Diagnostics     = table();
    
    month_in_year = Target_AllDate_date(floor(Target_AllDate_date / 10000) == year_now);

    for d = 1:length(month_in_year)
        
        tic;
    
        % Basic target information
        Target_Date = month_in_year(d);

        idx_target = Target_AllDate.date == Target_Date;

        Target_TTM    = Target_AllDate.TTM(idx_target);
        Target_ExDate = Target_AllDate.exdate(idx_target);

        % Diagnostics row
        Diag = struct();
        Diag.Target_Date          = Target_Date;
        Diag.Target_ExDate        = NaN;
        Diag.Target_TTM_Calendar  = NaN;
        Diag.Final_TTM_AfterAM    = NaN;
        Diag.Option_File          = "";
        Diag.N_Loaded             = NaN;
        Diag.N_SecID              = NaN;
        Diag.N_TargetDate         = NaN;
        Diag.N_TargetDateExDate   = NaN;
        Diag.N_AfterBidFilter     = NaN;
        Diag.N_AfterAskBidFilter  = NaN;
        Diag.N_AfterNoArb         = NaN;
        Diag.N_AfterOTM           = NaN;
        Diag.N_AfterDropNaNIV     = NaN;
        Diag.N_FinalObs           = NaN;
        Diag.N_UniqueK            = NaN;
        Diag.S0                   = NaN;
        Diag.S0_ADJ               = NaN;
        Diag.RF                   = NaN;
        Diag.DY                   = NaN;
        Diag.DY_Date_Used         = NaN;
        Diag.Mass_K               = NaN;
        Diag.Mass_R               = NaN;
        Diag.Neg_Ratio            = NaN;
        Diag.Min_RND              = NaN;
        Diag.Max_RND              = NaN;
        Diag.EQ_R                 = NaN;
        Diag.GEV_R_Success        = false;
        Diag.GEV_L_Success        = false;
        Diag.Status               = "Start";
        Diag.ElapsedSeconds       = NaN;

        if numel(Target_TTM) ~= 1
            warning('Target_Date %d has %d matching rows in Target_AllDate.', ...
                    Target_Date, numel(Target_TTM));

            Diag.Status = "Invalid target date mapping";
            Diag.ElapsedSeconds = toc;
            Table_Diagnostics = [Table_Diagnostics; struct2table(Diag, 'AsArray', true)];
            continue
        end

        Diag.Target_ExDate       = Target_ExDate;
        Diag.Target_TTM_Calendar = Target_TTM;

        disp(['Processing date: ', num2str(Target_Date)]);

    
        %*********************************************************************%
        % Step 1: Load Data
        % [1.  SecID              | 2.  Date (YYYYMMDD)    | 3.  TTM (Days)    | 4.  CPflag | 5.  K  | 6. S  
        % [7.  Option Price (Bid) | 8.  Option Price (Ask) | 9.  Open Interest | 10. Volume | 11. IV | 
        % [12. Delta              | 13. Gamma              | 14. Theta         | 15. Vega   ]        
        FileName = ['OP', ...
                    num2str(fix(Target_Date / 10000)), '_', ...
                    num2str(fix(rem(Target_Date, 10000) / 100)), '.txt'];
        FullFileName = fullfile(Path_Data_inc, FileName);
        Diag.Option_File = string(FileName);

        if exist(FullFileName, 'file') ~= 2
            warning('Option file does not exist: %s', FullFileName);

            Diag.Status = "Option file not found";
            Diag.ElapsedSeconds = toc;
            Table_Diagnostics = [Table_Diagnostics; struct2table(Diag, 'AsArray', true)];
            continue
        end

        Data = load(FullFileName);
        Diag.N_Loaded = size(Data, 1);

        clear FileName FullFileName
            
        % Index of Data
        Index_ID       = 1;
        Index_Date     = 2;
        Index_TTM      = 3;
        Index_CPFlag   = 4;
        Index_K        = 5;
        Index_S        = 6;
        Index_OP_Bid   = 7;
        Index_OP_Ask   = 8;
        Index_OI       = 9;
        Index_V        = 10;
        Index_IV       = 11;
        Index_Delta    = 12;
        Index_Gamma    = 13;
        Index_Theta    = 14;
        Index_Vega     = 15;

        % Constructed columns
        Index_RF       = Index_Vega + 1;
        Index_DY       = Index_Vega + 2;
        Index_S_ADJ    = Index_Vega + 3;


        %*********************************************************************%
        %  Step 2: Select SPX monthly option by Target_Date + Target_ExDate
        % 邏輯：raw expiration = quote date + raw TTM
        % 若 raw expiration 是 Saturday，先改成 Friday
        % 若該月 third Friday 是 Good Friday 等假日，target exdate 會是 Thursday
        % 再根據 Target_ExDate 修正成真正的 settlement date。

        Date_Num = datenum(num2str(Data(:, Index_Date)), 'yyyymmdd');

        Exp_Raw_Num     = Date_Num + Data(:, Index_TTM);
        Exp_Raw_Weekday = weekday(Exp_Raw_Num);

        % Saturday to Friday
        Exp_PreAM_Num = Exp_Raw_Num;
        Exp_PreAM_Num(Exp_Raw_Weekday == 7) = Exp_PreAM_Num(Exp_Raw_Weekday == 7) - 1;

        Target_Exp_Num = datenum(num2str(Target_ExDate), 'yyyymmdd');

        Exp_Corrected_Num = Exp_PreAM_Num;

        % Good Friday / holiday shift:
        % 若 target exdate 是週四，而 raw monthly contract 先修到週五，
        % 代表這個 Friday 是假日，應再往前修到 target exdate。
        idx_holiday_shifted_monthly = ...
            (weekday(Target_Exp_Num) ~= 6) & ...
            (Exp_Raw_Weekday == 7) & ...
            (Exp_PreAM_Num == Target_Exp_Num + 1);

        Exp_Corrected_Num(idx_holiday_shifted_monthly) = Target_Exp_Num;

        % Count before selection
        idx_id   = Data(:, Index_ID) == 108105;
        idx_date = idx_id & (Data(:, Index_Date) == Target_Date);

        Diag.N_SecID      = sum(idx_id);
        Diag.N_TargetDate = sum(idx_date);

        % Primary filter:
        % 1. SecID = 108105
        % 2. quote date = Target_Date
        % 3. raw expiration 是 Saturday：較能辨識 traditional monthly SPX
        % 4. corrected expiration = Target_ExDate
        idx_select = ...
            idx_id & ...
            (Data(:, Index_Date) == Target_Date) & ...
            (Exp_Raw_Weekday == 7) & ...
            (Exp_Corrected_Num == Target_Exp_Num);

        % Fallback:
        % 若某些檔案把 monthly 直接記 Friday / Thursday，不是 Saturday，
        % 放寬 raw Saturday 條件。
        if ~any(idx_select)

            warning('No Saturday-expiration monthly data found for %d. Relaxing monthly filter.', ...
                    Target_Date);

            idx_select = ...
                idx_id & ...
                (Data(:, Index_Date) == Target_Date) & ...
                (Exp_Corrected_Num == Target_Exp_Num);
        end

        Diag.N_TargetDateExDate = sum(idx_select);

        if ~any(idx_select)

            warning('No option data found for Target_Date = %d, Target_ExDate = %d.', ...
                    Target_Date, Target_ExDate);

            Diag.Status = "No data after Target_Date + Target_ExDate filter";
            Diag.ElapsedSeconds = toc;
            Table_Diagnostics = [Table_Diagnostics; struct2table(Diag, 'AsArray', true)];

            clear Data Date_Num Exp_Raw_Num Exp_Raw_Weekday Exp_PreAM_Num
            clear Target_Exp_Num Exp_Corrected_Num idx_holiday_shifted_monthly
            clear idx_id idx_date idx_select

            continue
        end

        Data = Data(idx_select, :);
        Date_Num_Selected = Date_Num(idx_select);
        Exp_Corrected_Num_Selected = Exp_Corrected_Num(idx_select);

        % Update TTM to corrected expiration before AM settlement
        Data(:, Index_TTM) = Exp_Corrected_Num_Selected - Date_Num_Selected;

        clear Date_Num Exp_Raw_Num Exp_Raw_Weekday Exp_PreAM_Num
        clear Target_Exp_Num Exp_Corrected_Num idx_holiday_shifted_monthly
        clear idx_id idx_date idx_select Date_Num_Selected Exp_Corrected_Num_Selected

    
        %*********************************************************************%
        % Step 3: AM Settlement Correction
        Data(:, Index_TTM) = Data(:, Index_TTM) - 1;
        Diag.Final_TTM_AfterAM = Data(1, Index_TTM);


        %*********************************************************************%
        % Step 4: Risk-Free Rate

        Data(:, Index_RF) = Function__RF_TTM(Data_RF, ...
                                             Data(:, Index_Date), ...
                                             Data(:, Index_TTM));

        if any(~isfinite(Data(:, Index_RF)))
            warning('Missing RF after interpolation for Target_Date = %d. Skip.', ...
                    Target_Date);

            Diag.Status = "Missing RF";
            Diag.ElapsedSeconds = toc;
            Table_Diagnostics = [Table_Diagnostics; struct2table(Diag, 'AsArray', true)];
            clear Data
            continue
        end


        %*********************************************************************%
        % Step 5: Dividend Yield

        Target_DateNum = datenum(num2str(Target_Date), 'yyyymmdd');
        idx_dy = (Data_DY(:, 1) == Target_Date);

        if ~any(idx_dy)

            % 若 target date 沒有 DY，往前找最多 10 天
            idx_dy_window = (Data_DY_DateNum <= Target_DateNum) & ...
                            (Data_DY_DateNum >= Target_DateNum - 10);

            if any(idx_dy_window)
                DY_Date_Used = max(Data_DY(idx_dy_window, 1));
                idx_dy = Data_DY(:, 1) == DY_Date_Used;
            else
                warning('No dividend yield found for Target_Date = %d.', Target_Date);

                Diag.Status = "Missing dividend yield";
                Diag.ElapsedSeconds = toc;
                Table_Diagnostics = [Table_Diagnostics; struct2table(Diag, 'AsArray', true)];
                clear Data Target_DateNum idx_dy idx_dy_window
                continue
            end

        else
            DY_Date_Used = Target_Date;
        end

        DY_Target = mean(Data_DY(idx_dy, end), 'omitnan');

        if ~isfinite(DY_Target)
            warning('Dividend yield is NaN for Target_Date = %d.', Target_Date);

            Diag.Status = "Dividend yield is NaN";
            Diag.ElapsedSeconds = toc;
            Table_Diagnostics = [Table_Diagnostics; struct2table(Diag, 'AsArray', true)];
            clear Data Target_DateNum idx_dy DY_Target
            continue
        end

        Data(:, Index_DY) = DY_Target;

        Diag.DY = DY_Target;
        Diag.DY_Date_Used = DY_Date_Used;

        clear Target_DateNum idx_dy DY_Target DY_Date_Used


        %*********************************************************************%
        % Step 6: Dividend-adjusted Stock Price for Black-Scholes

        S0  = Data(:, Index_S);
        TTM = Data(:, Index_TTM) / 365;
        DY  = Data(:, Index_DY);

        Data(:, Index_S_ADJ) = exp(-DY .* TTM) .* S0;

        clear S0 TTM DY

        
        %*********************************************************************%
        % Step 7: Data Filtering

        % 7.1 Bid Price > 3/8
        idx = Data(:, Index_OP_Bid) > (3 / 8);
        Data = Data(idx, :);
        Diag.N_AfterBidFilter = size(Data, 1);
        clear idx

        if isempty(Data)
            Diag.Status = "No data after bid filter";
            Diag.ElapsedSeconds = toc;
            Table_Diagnostics = [Table_Diagnostics; struct2table(Diag, 'AsArray', true)];
            continue
        end

        % 7.2 Ask Price > Bid Price > 0
        idx = (Data(:, Index_OP_Ask) > Data(:, Index_OP_Bid)) & ...
              (Data(:, Index_OP_Bid) > 0);

        Data = Data(idx, :);
        Diag.N_AfterAskBidFilter = size(Data, 1);
        clear idx

        if isempty(Data)
            Diag.Status = "No data after ask-bid filter";
            Diag.ElapsedSeconds = toc;
            Table_Diagnostics = [Table_Diagnostics; struct2table(Diag, 'AsArray', true)];
            continue
        end

        % 7.3 Standard no-arbitrage conditions
        K      = Data(:, Index_K);
        OP     = 0.5 * (Data(:, Index_OP_Bid) + Data(:, Index_OP_Ask));
        S0_ADJ = Data(:, Index_S_ADJ);
        TTM    = Data(:, Index_TTM) / 365;
        RF     = Data(:, Index_RF);

        idx_C = (Data(:, Index_CPFlag) == 1) & ...
                (OP <= S0_ADJ) & ...
                (OP >= S0_ADJ - exp(-RF .* TTM) .* K);

        idx_P = (Data(:, Index_CPFlag) == 2) & ...
                (OP <= exp(-RF .* TTM) .* K) & ...
                (OP >= exp(-RF .* TTM) .* K - S0_ADJ);

        idx = idx_C | idx_P;

        Data = Data(idx, :);
        Diag.N_AfterNoArb = size(Data, 1);

        clear K OP S0_ADJ TTM RF idx_C idx_P idx

        if isempty(Data)
            Diag.Status = "No data after no-arbitrage filter";
            Diag.ElapsedSeconds = toc;
            Table_Diagnostics = [Table_Diagnostics; struct2table(Diag, 'AsArray', true)];
            continue
        end

    
        %*********************************************************************%
        % Step 8: OTM and ATM blending

        K  = Data(:, Index_K);
        S0 = Data(:, Index_S);

        BP = 20;  % Figlewski (2010)

        idx = (Data(:, Index_CPFlag) == 1) & ...
              (K >= S0 - BP) & ...
              (K <= S0 + BP);

        Data(idx, Index_CPFlag) = 31;
        clear idx

        idx = (Data(:, Index_CPFlag) == 2) & ...
              (K >= S0 - BP) & ...
              (K <= S0 + BP);

        Data(idx, Index_CPFlag) = 32;
        clear idx BP

        idx = ((Data(:, Index_CPFlag) == 1) & (K >= S0)) | ...
              ((Data(:, Index_CPFlag) == 2) & (K <= S0)) | ...
              (fix(Data(:, Index_CPFlag) / 10) == 3);

        Data = Data(idx, :);
        Diag.N_AfterOTM = size(Data, 1);

        clear K S0 idx

        if isempty(Data)
            Diag.Status = "No data after OTM filter";
            Diag.ElapsedSeconds = toc;
            Table_Diagnostics = [Table_Diagnostics; struct2table(Diag, 'AsArray', true)];
            continue
        end

    
        %*********************************************************************%
        % Step 9: Step 9: Drop NaN IV
    
        idx = isfinite(Data(:, Index_IV));
        Data = Data(idx, :);

        Diag.N_AfterDropNaNIV = size(Data, 1);

        clear idx

        if isempty(Data)
            Diag.Status = "No data after dropping NaN IV";
            Diag.ElapsedSeconds = toc;
            Table_Diagnostics = [Table_Diagnostics; struct2table(Diag, 'AsArray', true)];
            continue
        end

        
        %*********************************************************************%
        % Step 10: Blend IV around ATM

        idx_bp = (Data(:, Index_CPFlag) == 31) | ...
                 (Data(:, Index_CPFlag) == 32);

        [AllK, ~, Index_AllK] = unique(Data(idx_bp, Index_K));

        if ~isempty(AllK)

            if length(AllK) > 1

                idx_bp_all = find(idx_bp);

                for i = 1:length(AllK)

                    idx_this = idx_bp_all(Index_AllK == i);
                    Data_BP = Data(idx_this, :);
                    Data_BP = sortrows(Data_BP, -Index_CPFlag);            % Put to Call

                    if size(Data_BP, 1) > 1

                        W = (max(AllK) - AllK(i)) / (max(AllK) - min(AllK));
                        Data(idx_this, Index_IV) = [W, 1 - W] * Data_BP(:, Index_IV);

                    else

                        Data(idx_this, Index_CPFlag) = 3;

                    end

                    clear idx_this Data_BP W
                end

                clear idx_bp_all i

            else

                idx_bp_all = find(idx_bp);

                if length(idx_bp_all) > 1
                    Data(idx_bp_all, Index_IV) = mean(Data(idx_bp_all, Index_IV), 'omitnan');
                else
                    Data(idx_bp_all, Index_CPFlag) = 3;
                end

                clear idx_bp_all
            end
        end

        clear idx_bp AllK Index_AllK

        Data(Data(:, Index_CPFlag) == 31, Index_CPFlag) = 3;
        Data(Data(:, Index_CPFlag) == 32, :) = [];

        Data = sortrows(Data, Index_K);

    
        %*********************************************************************%
        % Step 11: Check data sufficiency before spline

        Num_Obs     = size(Data, 1);
        Num_UniqueK = numel(unique(Data(:, Index_K)));

        Diag.N_FinalObs = Num_Obs;
        Diag.N_UniqueK  = Num_UniqueK;

        if Num_Obs < 10 || Num_UniqueK < 6
            warning('Too few option observations for Target_Date = %d. Obs = %d, Unique K = %d. Skip.', ...
                    Target_Date, Num_Obs, Num_UniqueK);

            Diag.Status = "Too few observations";
            Diag.ElapsedSeconds = toc;
            Table_Diagnostics = [Table_Diagnostics; struct2table(Diag, 'AsArray', true)];
            clear Data
            continue
        end

    
        %*********************************************************************%
        % Step 12: Determine RND range and smooth IV

        S0      = Data(1, Index_S);
        S0_ADJ0 = Data(1, Index_S_ADJ);
        TTM0    = Data(1, Index_TTM) / 365;
        RF0     = Data(1, Index_RF);
        DY0     = Data(1, Index_DY);

        Diag.S0     = S0;
        Diag.S0_ADJ = S0_ADJ0;
        Diag.RF     = RF0;
        Diag.DY     = DY0;

        K_Low  = (0.3 / 100) * S0;
        K_High = 3 * S0;

        Smooth_AllK = linspace(K_Low, K_High, Num_Grid)';

        clear K_Low K_High

        % Use unique K and average IV to avoid spline instability
        [K_Unique, ~, idx_K_Group] = unique(Data(:, Index_K));
        IV_Unique = accumarray(idx_K_Group, Data(:, Index_IV), [], @mean);

        try
            LSS = spap2(1, 4, K_Unique, IV_Unique);
            LSS = spap2(newknt(LSS), 4, K_Unique, IV_Unique);
        catch ME
            warning('Spline fitting failed for Target_Date = %d. Message: %s', ...
                    Target_Date, ME.message);

            Diag.Status = "Spline fitting failed";
            Diag.ElapsedSeconds = toc;
            Table_Diagnostics = [Table_Diagnostics; struct2table(Diag, 'AsArray', true)];

            clear Data Smooth_AllK K_Unique IV_Unique idx_K_Group
            continue
        end

        Smooth_K = Smooth_AllK((Smooth_AllK >= min(K_Unique)) & ...
                               (Smooth_AllK <= max(K_Unique)));

        Smooth_IV = fnval(LSS, Smooth_K);
        Smooth_IV = Smooth_IV(:);

        clear LSS K_Unique IV_Unique idx_K_Group

        % Check valid smoothed IV
        idx_valid_iv = isfinite(Smooth_IV) & (Smooth_IV > 0);

        if sum(idx_valid_iv) < 100
            warning('Too few valid smoothed IV values for Target_Date = %d. Skip.', ...
                    Target_Date);

            Diag.Status = "Too few valid smoothed IV";
            Diag.ElapsedSeconds = toc;
            Table_Diagnostics = [Table_Diagnostics; struct2table(Diag, 'AsArray', true)];

            clear Data Smooth_AllK Smooth_K Smooth_IV idx_valid_iv
            continue
        end

        Smooth_K  = Smooth_K(idx_valid_iv);
        Smooth_IV = Smooth_IV(idx_valid_iv);

        clear idx_valid_iv

    
        %*********************************************************************%
        % Step 13: Smoothed option price and empirical RND

        S0_ADJ_vec = S0_ADJ0 * ones(size(Smooth_K));
        TTM_vec    = TTM0    * ones(size(Smooth_K));
        RF_vec     = RF0     * ones(size(Smooth_K));

        Smooth_OP = blsprice(S0_ADJ_vec, Smooth_K, RF_vec, TTM_vec, Smooth_IV, []);

        if any(~isfinite(Smooth_OP))
            warning('Smoothed option price contains NaN/Inf for Target_Date = %d. Skip.', ...
                    Target_Date);

            Diag.Status = "Invalid smoothed option price";
            Diag.ElapsedSeconds = toc;
            Table_Diagnostics = [Table_Diagnostics; struct2table(Diag, 'AsArray', true)];

            clear Data Smooth_AllK Smooth_K Smooth_IV Smooth_OP
            clear S0_ADJ_vec TTM_vec RF_vec
            continue
        end

        % Empirical PDF and CDF on K_2 to K_{N-1}
        Smooth_EMP_PDF = exp(RF0 * TTM0) * ...
                         (Smooth_OP(3:end) - 2 * Smooth_OP(2:end-1) + Smooth_OP(1:end-2)) ./ ...
                         ((Smooth_K(3:end) - Smooth_K(2:end-1)) .^ 2);

        Smooth_EMP_CDF = exp(RF0 * TTM0) * ...
                         ((Smooth_OP(3:end) - Smooth_OP(1:end-2)) ./ ...
                         (Smooth_K(3:end) - Smooth_K(1:end-2))) + 1;

        Smooth_K  = Smooth_K(2:end-1);
        Smooth_OP = Smooth_OP(2:end-1);
        Smooth_IV = Smooth_IV(2:end-1);

        Smooth_EMP_PDF = Smooth_EMP_PDF(:);
        Smooth_EMP_CDF = Smooth_EMP_CDF(:);

        idx_emp_valid = isfinite(Smooth_EMP_PDF) & ...
                        isfinite(Smooth_EMP_CDF) & ...
                        (Smooth_EMP_PDF > 0) & ...
                        (Smooth_EMP_CDF > 0) & ...
                        (Smooth_EMP_CDF < 1);

        if sum(idx_emp_valid) < 10
            warning('Too few valid empirical PDF/CDF points for Target_Date = %d. Skip.', ...
                    Target_Date);

            Diag.Status = "Too few valid empirical PDF/CDF";
            Diag.ElapsedSeconds = toc;
            Table_Diagnostics = [Table_Diagnostics; struct2table(Diag, 'AsArray', true)];

            clear Data Smooth_AllK Smooth_K Smooth_IV Smooth_OP
            clear Smooth_EMP_PDF Smooth_EMP_CDF idx_emp_valid
            clear S0_ADJ_vec TTM_vec RF_vec
            continue
        end

        Smooth_K_Valid       = Smooth_K(idx_emp_valid);
        Smooth_EMP_PDF_Valid = Smooth_EMP_PDF(idx_emp_valid);
        Smooth_EMP_CDF_Valid = Smooth_EMP_CDF(idx_emp_valid);

        clear idx_emp_valid S0_ADJ_vec TTM_vec RF_vec

    
        %*********************************************************************%
        % Step 14: Right-tail GEV connection

        Smooth_GEV_R_PDF = nan(size(Smooth_AllK));
        Smooth_GEV_R_CDF = nan(size(Smooth_AllK));
        Parameters_GEV_R = nan(1, 3);
        FitError_GEV_R   = nan(1, 3);
        flag_GEV_R       = false;

        BP_R1 = 0.98;

        idx = find((Smooth_EMP_CDF_Valid - BP_R1) >= 0, 1, 'first');

        if ~isempty(idx)
            BP_R1      = Smooth_EMP_CDF_Valid(idx);
            K_R1       = Smooth_K_Valid(idx);
            EMP_PDF_R1 = Smooth_EMP_PDF_Valid(idx);
        else
            BP_R1      = Smooth_EMP_CDF_Valid(end);
            K_R1       = Smooth_K_Valid(end);
            EMP_PDF_R1 = Smooth_EMP_PDF_Valid(end);
        end

        BP_R0 = BP_R1 - 0.03;

        idx = find((Smooth_EMP_CDF_Valid - BP_R0) >= 0, 1, 'first');

        if isempty(idx)
            warning('Cannot find right-tail BP_R0 for Target_Date = %d. Skip.', Target_Date);

            Diag.Status = "Cannot find right-tail BP_R0";
            Diag.ElapsedSeconds = toc;
            Table_Diagnostics = [Table_Diagnostics; struct2table(Diag, 'AsArray', true)];

            clear Data Smooth_AllK Smooth_K Smooth_IV Smooth_OP
            clear Smooth_EMP_PDF Smooth_EMP_CDF Smooth_K_Valid Smooth_EMP_PDF_Valid Smooth_EMP_CDF_Valid
            continue
        end

        BP_R0      = Smooth_EMP_CDF_Valid(idx);
        K_R0       = Smooth_K_Valid(idx);
        EMP_PDF_R0 = Smooth_EMP_PDF_Valid(idx);
        EMP_CDF_R0 = Smooth_EMP_CDF_Valid(idx);

        clear idx

        try
            [mu, sigma, k] = Function__Fit_GEV_PDF(K_R0, K_R1, ...
                                                   EMP_CDF_R0, EMP_PDF_R0, EMP_PDF_R1);

            Smooth_GEV_R_PDF = gevpdf(Smooth_AllK, k, sigma, mu);
            Smooth_GEV_R_CDF = gevcdf(Smooth_AllK, k, sigma, mu);

            Parameters_GEV_R = [mu, sigma, k];

            FitError_GEV_R = Function__CheckError_GEV_PDFCDF(mu, sigma, k, ...
                                                             K_R0, K_R1, ...
                                                             EMP_CDF_R0, EMP_PDF_R0, EMP_PDF_R1);

            flag_GEV_R = true;

        catch ME
            warning('Right-tail GEV fitting failed for Target_Date = %d. Message: %s', ...
                    Target_Date, ME.message);
        end

        Diag.GEV_R_Success = flag_GEV_R;

        clear mu sigma k EMP_CDF_R0 EMP_PDF_R0 EMP_PDF_R1

    
        %*********************************************************************%
        % Step 15: Left-tail GEV connection

        Smooth_GEV_L_PDF = nan(size(Smooth_AllK));
        Smooth_GEV_L_CDF = nan(size(Smooth_AllK));
        Parameters_GEV_L = nan(1, 3);
        FitError_GEV_L   = nan(1, 3);
        flag_GEV_L       = false;

        BP_L1 = 0.02;

        idx = find((BP_L1 - Smooth_EMP_CDF_Valid) >= 0, 1, 'last');

        if ~isempty(idx)
            BP_L1      = Smooth_EMP_CDF_Valid(idx);
            K_L1       = Smooth_K_Valid(idx);
            EMP_PDF_L1 = Smooth_EMP_PDF_Valid(idx);
        else
            BP_L1      = Smooth_EMP_CDF_Valid(1);
            K_L1       = Smooth_K_Valid(1);
            EMP_PDF_L1 = Smooth_EMP_PDF_Valid(1);
        end

        BP_L0 = BP_L1 + 0.03;

        idx = find((BP_L0 - Smooth_EMP_CDF_Valid) >= 0, 1, 'last');

        if isempty(idx)
            warning('Cannot find left-tail BP_L0 for Target_Date = %d. Skip.', Target_Date);

            Diag.Status = "Cannot find left-tail BP_L0";
            Diag.ElapsedSeconds = toc;
            Table_Diagnostics = [Table_Diagnostics; struct2table(Diag, 'AsArray', true)];

            clear Data Smooth_AllK Smooth_K Smooth_IV Smooth_OP
            clear Smooth_EMP_PDF Smooth_EMP_CDF Smooth_K_Valid Smooth_EMP_PDF_Valid Smooth_EMP_CDF_Valid
            continue
        end

        BP_L0      = Smooth_EMP_CDF_Valid(idx);
        K_L0       = Smooth_K_Valid(idx);
        EMP_PDF_L0 = Smooth_EMP_PDF_Valid(idx);
        EMP_CDF_L0 = Smooth_EMP_CDF_Valid(idx);

        clear idx

        try
            [mu, sigma, k] = Function__Fit_GEV_PDF(-K_L0, -K_L1, ...
                                                   1 - EMP_CDF_L0, EMP_PDF_L0, EMP_PDF_L1);

            Smooth_GEV_L_PDF = gevpdf(-Smooth_AllK, k, sigma, mu);
            Smooth_GEV_L_CDF = 1 - gevcdf(-Smooth_AllK, k, sigma, mu);

            Parameters_GEV_L = [-mu, sigma, k];

            FitError_GEV_L = Function__CheckError_GEV_PDFCDF(mu, sigma, k, ...
                                                             -K_L0, -K_L1, ...
                                                             1 - EMP_CDF_L0, EMP_PDF_L0, EMP_PDF_L1);

            flag_GEV_L = true;

        catch ME
            warning('Left-tail GEV fitting failed for Target_Date = %d. Message: %s', ...
                    Target_Date, ME.message);
        end

        Diag.GEV_L_Success = flag_GEV_L;

        clear mu sigma k EMP_CDF_L0 EMP_PDF_L0 EMP_PDF_L1

    
        %*********************************************************************%
        % Step 16: Combine full RND

        if ~flag_GEV_R || ~flag_GEV_L
            warning('GEV tail fitting failed for Target_Date = %d. Skip this month.', ...
                    Target_Date);

            Diag.Status = "GEV fitting failed";
            Diag.ElapsedSeconds = toc;
            Table_Diagnostics = [Table_Diagnostics; struct2table(Diag, 'AsArray', true)];

            clear Data Smooth_AllK Smooth_K Smooth_IV Smooth_OP
            clear Smooth_EMP_PDF Smooth_EMP_CDF Smooth_K_Valid Smooth_EMP_PDF_Valid Smooth_EMP_CDF_Valid
            clear Smooth_GEV_R_PDF Smooth_GEV_R_CDF Smooth_GEV_L_PDF Smooth_GEV_L_CDF
            continue
        end

        Smooth_AllK_RND = nan(size(Smooth_AllK));

        % Empirical center: use interpolation to avoid length mismatch
        idx_emp = (Smooth_AllK >= K_L0) & ...
                  (Smooth_AllK <= K_R0);

        Smooth_AllK_RND(idx_emp) = interp1(Smooth_K, Smooth_EMP_PDF, ...
                                           Smooth_AllK(idx_emp), ...
                                           'linear');

        % Right tail
        idx_gev_r = Smooth_AllK >= K_R0;
        Smooth_AllK_RND(idx_gev_r) = Smooth_GEV_R_PDF(idx_gev_r);

        % Left tail
        idx_gev_l = Smooth_AllK <= K_L0;
        Smooth_AllK_RND(idx_gev_l) = Smooth_GEV_L_PDF(idx_gev_l);

        clear idx_emp idx_gev_r idx_gev_l

    
        %*********************************************************************%
        % Step 17: RND cleaning and normalization

        Smooth_AllK_RND = Smooth_AllK_RND(:);
        Smooth_AllK     = Smooth_AllK(:);

        if any(~isfinite(Smooth_AllK_RND))
            warning('RND contains NaN or Inf for Target_Date = %d. Skip this month.', ...
                    Target_Date);

            Diag.Status = "RND contains NaN or Inf";
            Diag.ElapsedSeconds = toc;
            Table_Diagnostics = [Table_Diagnostics; struct2table(Diag, 'AsArray', true)];

            clear Data Smooth_AllK Smooth_AllK_RND
            continue
        end

        % 計算負密度面積與正密度面積
        neg_mass = trapz(Smooth_AllK, max(-Smooth_AllK_RND, 0));
        pos_mass = trapz(Smooth_AllK, max( Smooth_AllK_RND, 0));

        if pos_mass <= 0 || ~isfinite(pos_mass)
            warning('RND has non-positive positive mass for Target_Date = %d. Skip.', ...
                    Target_Date);

            Diag.Status = "RND has non-positive mass";
            Diag.ElapsedSeconds = toc;
            Table_Diagnostics = [Table_Diagnostics; struct2table(Diag, 'AsArray', true)];

            clear Data Smooth_AllK Smooth_AllK_RND
            continue
        end

        neg_ratio = neg_mass / pos_mass;

        if neg_ratio > 0.01
            warning('Large negative RND mass for Target_Date = %d. Negative mass ratio = %.4f.', ...
                    Target_Date, neg_ratio);
        end

        % Truncate negative density for numerical stability
        Smooth_AllK_RND(Smooth_AllK_RND < 0) = 0;

        % 正規化：K-space density 積分成 1
        mass_K_before = trapz(Smooth_AllK, Smooth_AllK_RND);

        if ~isfinite(mass_K_before) || mass_K_before <= 0
            warning('Invalid RND mass for Target_Date = %d. Skip.', Target_Date);

            Diag.Status = "Invalid RND mass";
            Diag.ElapsedSeconds = toc;
            Table_Diagnostics = [Table_Diagnostics; struct2table(Diag, 'AsArray', true)];

            clear Data Smooth_AllK Smooth_AllK_RND
            continue
        end

        Smooth_AllK_RND = Smooth_AllK_RND / mass_K_before;

        mass_K_after = trapz(Smooth_AllK, Smooth_AllK_RND);

        if abs(mass_K_after - 1) > 1e-6
            warning('K-space RND normalization failed for Target_Date = %d. Mass_K = %.10f.', ...
                    Target_Date, mass_K_after);
        end

    
        %*********************************************************************%
        % Step 18: Transform from K-space to R-space

        % 採用 Linn / Schreindorfer 的 ex-dividend return 做法：
        % R = K / S0
        % f_R(R) = S0 * f_K(K)
        % 因此這裡不用 S0_ADJ
        % E^Q[R] 小於 Rf 是可以接受的，不在這裡強制修正。

        Smooth_AllR     = Smooth_AllK / S0;
        Smooth_AllR_RND = S0 * Smooth_AllK_RND;

        mass_R = trapz(Smooth_AllR, Smooth_AllR_RND);

        if abs(mass_R - 1) > 1e-6
            warning('R-space RND does not integrate to one for Target_Date = %d. Mass_R = %.10f.', ...
                    Target_Date, mass_R);
        end

        EQ_R = trapz(Smooth_AllR, Smooth_AllR_RND .* Smooth_AllR);

    
        %*********************************************************************%
        % Step 19: Output tables

        columnName = num2str(Target_Date);

        Table_Smooth_AllK.(columnName)     = Smooth_AllK;
        Table_Smooth_AllR.(columnName)     = Smooth_AllR;
        Table_Smooth_AllR_RND.(columnName) = Smooth_AllR_RND;

        Diag.Mass_K    = mass_K_after;
        Diag.Mass_R    = mass_R;
        Diag.Neg_Ratio = neg_ratio;
        Diag.Min_RND   = min(Smooth_AllK_RND);
        Diag.Max_RND   = max(Smooth_AllK_RND);
        Diag.EQ_R      = EQ_R;
        Diag.Status    = "Success";
        Diag.ElapsedSeconds = toc;

        Table_Diagnostics = [Table_Diagnostics; struct2table(Diag, 'AsArray', true)];

        disp(['     Spend Time: ', num2str(Diag.ElapsedSeconds), ' Seconds']);

        clear Data
        clear Smooth_AllK Smooth_AllR Smooth_AllK_RND Smooth_AllR_RND
        clear Smooth_K Smooth_IV Smooth_OP Smooth_EMP_PDF Smooth_EMP_CDF
        clear Smooth_K_Valid Smooth_EMP_PDF_Valid Smooth_EMP_CDF_Valid
        clear Smooth_GEV_R_PDF Smooth_GEV_R_CDF Parameters_GEV_R FitError_GEV_R
        clear Smooth_GEV_L_PDF Smooth_GEV_L_CDF Parameters_GEV_L FitError_GEV_L
        clear BP_R0 BP_R1 K_R0 K_R1 BP_L0 BP_L1 K_L0 K_L1
        clear S0 S0_ADJ0 TTM0 RF0 DY0
        clear mass_K_before mass_K_after mass_R neg_mass pos_mass neg_ratio EQ_R
        clear flag_GEV_R flag_GEV_L columnName Diag

    end


    % =====================================================================
    %  Save yearly output

    FileName_MAT = fullfile(Path_Output, ...
        ['TTM_', num2str(Target_AllTTM), '_RND_Tables_', num2str(year_now), '.mat']);

    save(FileName_MAT, ...
         'Table_Smooth_AllK', ...
         'Table_Smooth_AllR', ...
         'Table_Smooth_AllR_RND', ...
         'Table_Diagnostics');

    FileName_Diag = fullfile(Path_Output, ...
        ['TTM_', num2str(Target_AllTTM), '_RND_Diagnostics_', num2str(year_now), '.csv']);

    writetable(Table_Diagnostics, FileName_Diag);

    Table_Diagnostics_All = [Table_Diagnostics_All; Table_Diagnostics];

    clear FileName_MAT FileName_Diag
    clear Table_Smooth_AllK Table_Smooth_AllR Table_Smooth_AllR_RND Table_Diagnostics
end


% ========================================================================
%  Save combined diagnostics

FileName_Diag_All = fullfile(Path_Output, ...
    ['TTM_', num2str(Target_AllTTM), '_RND_Diagnostics_All.csv']);

writetable(Table_Diagnostics_All, FileName_Diag_All);

disp('All done.');

clear Data_RF Data_DY Data_DY_DateNum
clear Table_Diagnostics_All FileName_Diag_All
clear Target_AllDate Target_AllDate_date Target_AllDate_exdate
clear Path_MainFolder Path_Data Path_Output Path_Data_inc Path_RND_Function
clear Target_AllTTM Num_Grid years y d year_now month_in_year