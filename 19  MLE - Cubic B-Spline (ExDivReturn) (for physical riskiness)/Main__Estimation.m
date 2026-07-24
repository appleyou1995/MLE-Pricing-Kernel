clear; clc;

%% Paths

Path_MainFolder = ['D:\Google\我的雲端硬碟\學術｜研究與論文\論文著作\' ...
    'MLE Pricing Kernel'];
Path_Data   = fullfile(Path_MainFolder, 'Code', '00  Output');
Path_RND    = fullfile(Path_MainFolder, 'Code', '01  Output');
Path_Code19 = fullfile(Path_MainFolder, 'Code', ...
    '19  MLE - Cubic B-Spline (ExDivReturn) (for physical riskiness)');
Path_Output = fullfile(Path_MainFolder, 'Code', '19  Output');

addpath(Path_Code19);
if ~exist(Path_Output, 'dir'), mkdir(Path_Output); end


%% Run configuration

TTM_List     = [30, 60, 90, 180];
Sample_Start = 19960101;
Sample_End   = 20241231;

Spline_Degree      = 3;     % Cubic B-spline.
Num_Basis_Function = 5;     % Cubic spline with exactly one internal knot.
Enforce_Decreasing = true;  % Preserve the decreasing-SDF restriction.
Expected_Grid_Size = 30000;

fprintf(['Cubic B-spline MLE without probability distortion\n' ...
    'TTM list: %s\nSample: %08d to %08d\n\n'], ...
    mat2str(TTM_List), Sample_Start, Sample_End);


%% Shared input

Risk_Free_Rate_All = readtable(fullfile( ...
    Path_Data, 'Risk_Free_GrossFactor_ByTargetTTM.csv'));

Theta_All   = table();
Kappa_All   = table();
Summary_All = table();
Manifest_All = table();


%% Estimate every TTM separately

for ttm_idx = 1:numel(TTM_List)
    Target_TTM = TTM_List(ttm_idx);
    fprintf('\n============================================================\n');
    fprintf('Starting TTM = %d\n', Target_TTM);
    fprintf('============================================================\n');

    Data = load_ttm_inputs( ...
        Path_Data, Path_RND, Risk_Free_Rate_All, Target_TTM, ...
        Sample_Start, Sample_End, Expected_Grid_Size);

    fprintf(['Aligned sample: %d months, %08d to %08d; ' ...
        'grid size = %d per month.\n'], ...
        numel(Data.Quote_Dates), Data.Quote_Dates(1), ...
        Data.Quote_Dates(end), Expected_Grid_Size);

    [theta_hat, log_lik, BIC, exitflag, optim_output, ...
        kappa_vec, SDF_Cell, ~, ~, Model_Info] = ...
        MLE_BSpline_estimation( ...
            Data.Smooth_AllR, Data.Smooth_AllR_RND, ...
            Data.Realized_Return, Data.Risk_Free_Rate, ...
            Num_Basis_Function, Spline_Degree, Enforce_Decreasing, ...
            Data.Global_Min_R, Data.Global_Max_R);

    [Theta_Table, Kappa_Table, Summary_Table, Manifest_Table] = ...
        export_ttm_results_csv( ...
            Path_Output, Target_TTM, Data, theta_hat, kappa_vec, ...
            SDF_Cell, log_lik, BIC, ...
            exitflag, optim_output, Model_Info);

    Theta_All    = [Theta_All; Theta_Table];       %#ok<AGROW>
    Kappa_All    = [Kappa_All, Kappa_Table];       %#ok<AGROW>
    Summary_All  = [Summary_All; Summary_Table];   %#ok<AGROW>
    Manifest_All = [Manifest_All; Manifest_Table]; %#ok<AGROW>

    clear Data theta_hat kappa_vec SDF_Cell
    clear Model_Info Theta_Table Kappa_Table Summary_Table Manifest_Table
end


%% Combined CSV files for Python handoff

writetable(Theta_All, fullfile(Path_Output, 'Theta_By_TTM.csv'));
writetable(Kappa_All, fullfile(Path_Output, 'Kappa_By_TTM.csv'));
writetable(Summary_All, fullfile(Path_Output, ...
    'CubicBSpline_EstimationSummary_AllTTM.csv'));
writetable(Manifest_All, fullfile(Path_Output, ...
    'CubicBSpline_SDFGrid_Manifest.csv'));

fprintf('\nAll four TTM estimations and CSV exports are complete.\n');
fprintf('Output folder: %s\n', Path_Output);
