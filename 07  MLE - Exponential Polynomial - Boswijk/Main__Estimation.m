clear; clc;

Path_MainFolder = 'D:\Google\我的雲端硬碟\學術｜研究與論文\論文著作\MLE Pricing Kernel';
Path_Data = 'D:\Google\我的雲端硬碟\學術｜研究與論文\論文著作\CDI Method';


%% Load the data

Target_TTM = 30;

% Load realized gross returns (R_{t+1})
Path_Data_01 = fullfile(Path_Data, 'Code', '01  輸出資料');
FileName = ['Realized_Return_TTM_', num2str(Target_TTM), '.csv'];
Realized_Return = readtable(fullfile(Path_Data_01, FileName));

% Load risk-free rate R_f^t
Path_Data_01_main = fullfile(Path_Data, 'Code', '01  原始資料處理');
FileName = 'Risk_Free_Rate.csv';
Risk_Free_Rate_All = readtable(fullfile(Path_Data_01_main, FileName));
Risk_Free_Rate = Risk_Free_Rate_All.rf_gross_29d;

% Load Q-measure PDF tables: R axis and corresponding f^*_t(R)
Path_Data_02 = fullfile(Path_Data, 'Code', '02  輸出資料');
Smooth_AllK = [];
Smooth_AllR = [];
Smooth_AllR_RND = [];

years_to_merge = 1996:2021;
for year = years_to_merge
    input_filename = fullfile(Path_Data_02, sprintf('TTM_%d_RND_Tables_%d.mat', Target_TTM, year));
    if exist(input_filename, 'file')
        data = load(input_filename);
        Smooth_AllK = [Smooth_AllK, data.Table_Smooth_AllK];
        Smooth_AllR = [Smooth_AllR, data.Table_Smooth_AllR];               % R_grid for interpolation
        Smooth_AllR_RND = [Smooth_AllR_RND, data.Table_Smooth_AllR_RND];   % f^*_t(R) on grid
    else
        warning('File %s does not exist.', input_filename);
    end
end

clear Path_Data_01 Path_Data_01_main Path_Data_02 Target_TTM
clear Risk_Free_Rate_All years_to_merge data FileName input_filename year


%% 

dates = intersect(Smooth_AllR.Properties.VariableNames, ...
                  Smooth_AllR_RND.Properties.VariableNames, 'stable');

for k = 1:numel(dates)
    d = dates{k};

    R = Smooth_AllR{1, d};        % 1x30000 row（不要用 (:)）
    F = Smooth_AllR_RND{1, d};    % 1x30000 row

    % 清理與正規化
    F(~isfinite(F)) = 0;
    F = max(F, 0);
    I = trapz(R, F);

    if isfinite(I) && I > 0
        Smooth_AllR_RND{1, d} = F / I;   % 寫回 1x30000 row
    else
        warning('RND:%s:BadIntegral','%s',d)
    end
end



%% 

function rpt = sanity_check_R_and_RND(Smooth_AllR, Smooth_AllR_RND, varargin)
% 扫描每個日期欄，檢查 R 軸與 RND 的常見問題，回傳表格報表。
% 參數：'TolIntegral' (default 1e-3), 'TolLogRRange' (default 50)

ip = inputParser;
ip.addParameter('TolIntegral', 1e-3);
ip.addParameter('TolLogRRange', 50);
ip.parse(varargin{:});
TolInt   = ip.Results.TolIntegral;
TolRange = ip.Results.TolLogRRange;

datesA = Smooth_AllR.Properties.VariableNames;
datesB = Smooth_AllR_RND.Properties.VariableNames;

% 只掃描兩邊皆存在的日期
dates = intersect(datesA, datesB, 'stable');

% 預先配置
n = numel(dates);
Date        = strings(n,1);
nR          = zeros(n,1);
nF          = zeros(n,1);
NonFiniteR  = false(n,1);
NonFiniteF  = false(n,1);
NonPosR     = false(n,1);
NotMonoR    = false(n,1);
DupR        = false(n,1);
NegFracF    = NaN(n,1);
IntF        = NaN(n,1);
IntDev      = NaN(n,1);
LogRRange   = NaN(n,1);

for k = 1:n
    d = dates{k};
    Date(k) = string(d);

    R = Smooth_AllR.(d);      R = R(:);     % 轉成 column
    F = Smooth_AllR_RND.(d);  F = F(:);

    nR(k) = numel(R);
    nF(k) = numel(F);

    % --- R 軸檢查 ---
    NonFiniteR(k) = any(~isfinite(R));
    NonPosR(k)    = any(R <= 0);
    DupR(k)       = (numel(unique(R)) < numel(R));
    NotMonoR(k)   = any(diff(R) <= 0);      % 要嚴格遞增才安全

    % logR 範圍（先夾正避免 log(<=0)）
    Rpos          = max(R, eps);
    logR          = log(Rpos);
    LogRRange(k)  = max(logR) - min(logR);

    % --- F 檢查 ---
    NonFiniteF(k) = any(~isfinite(F));
    if nR(k) >= 2 && nF(k) >= 2
        NegFracF(k) = mean(F < 0, 'omitnan');

        % 积分（把 F 的 NaN/負值先歸零以便觀察偏差）
        Fclean = F;
        Fclean(~isfinite(Fclean)) = 0;
        Fclean = max(Fclean, 0);
        IntF(k)   = trapz(R, Fclean);
        IntDev(k) = abs(IntF(k) - 1);
    end
end

% 彙整報表
rpt = table( ...
    Date, nR, nF, NonFiniteR, NonPosR, NotMonoR, DupR, LogRRange, ...
    NonFiniteF, NegFracF, IntF, IntDev);

% 加一個 Flag 欄位方便你直接挑異常行
rpt.Flag = (NonFiniteR | NonPosR | NotMonoR | DupR | ...
            NonFiniteF | (NegFracF > 0) | (IntDev > TolInt) | (LogRRange > TolRange));

% 友善列印：只顯示有問題的 Top 10
bad = rpt(rpt.Flag,:);
if ~isempty(bad)
    disp('⚠️ 發現異常的日期（前10筆）：');
    disp(bad(1:min(10,height(bad)), :));
else
    disp('✅ 沒發現異常。');
end
end


rpt = sanity_check_R_and_RND(Smooth_AllR, Smooth_AllR_RND, ...
       'TolIntegral', 1e-3, 'TolLogRRange', 50);

% 只看有問題的
bad = rpt(rpt.Flag,:);



%% (1) Estimation - with delta

Path_Output = fullfile(Path_MainFolder, 'Code', '07  Output - with delta');

% Add paths
Path_Code_07 = fullfile(Path_MainFolder, 'Code', ...
    '07  MLE - Exponential Polynomial - Boswijk');
addpath(Path_Code_07);

% Initialize result storage
max_L     = 3;
use_delta = true;

for L = 1:max_L
    fprintf('\n--- Estimating MLE with δ, max L = %d ---\n', L);

    [gamma_hat, log_lik, delta_vec, M_vec] = MLE_theta_estimation( ...
        Smooth_AllR, Smooth_AllR_RND, Realized_Return, Risk_Free_Rate, L, use_delta);

    OutputFile = fullfile(Path_Output, sprintf('MLE_gamma_max_L_%d.mat', L));
    save(OutputFile, 'gamma_hat', 'log_lik', 'delta_vec', 'M_vec');
end


%% (2) Estimation - without delta

Path_Output = fullfile(Path_MainFolder, 'Code', '07  Output - without delta');

% Add paths
Path_Code_07 = fullfile(Path_MainFolder, 'Code', ...
    '07  MLE - Exponential Polynomial - Boswijk');
addpath(Path_Code_07);

% Initialize result storage
max_L     = 3;
use_delta = false;

for L = 1:max_L
    fprintf('\n--- Estimating MLE without δ, max L = %d ---\n', L);

    [gamma_hat, log_lik, ~, M_vec] = MLE_theta_estimation( ...
        Smooth_AllR, Smooth_AllR_RND, Realized_Return, Risk_Free_Rate, L, use_delta);

    OutputFile = fullfile(Path_Output, sprintf('MLE_gamma_max_L_%d.mat', L));
    save(OutputFile, 'gamma_hat', 'log_lik', 'M_vec');
end