clear; clc;

%% Paths

Path_MainFolder = ['D:\Google\我的雲端硬碟\學術｜研究與論文\論文著作\' ...
    'MLE Pricing Kernel'];
Path_Data   = fullfile(Path_MainFolder, 'Code', '00  Hsieh');
Path_RND    = fullfile(Path_MainFolder, 'Code', '01  Output - Hsieh');
Path_Code19 = fullfile(Path_MainFolder, 'Code', ...
    '19  MLE - Cubic B-Spline (ExDivReturn) (for physical riskiness)');
Path_Output = fullfile(Path_MainFolder, 'Code', '19  Output - Hsieh');

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


%% Shared input and compact output containers

Risk_Free_Rate_All = readtable(fullfile( ...
    Path_Data, 'Risk_Free_GrossFactor_ByTargetTTM.csv'));

Theta_All      = table();
Summary_All    = table();
G_Shape_All    = table();
G_TimeScale_All = table();
R_Grid_All     = table();
Kappa_By_TTM   = cell(numel(TTM_List), 1);


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

    [Theta_Table, Kappa_Table, Summary_Table, R_Grid_Table, ...
        G_Shape_Table, G_TimeScale_Table] = export_ttm_results_csv( ...
            Target_TTM, Data, theta_hat, kappa_vec, SDF_Cell, ...
            log_lik, BIC, exitflag, optim_output, Model_Info);

    Theta_All       = [Theta_All; Theta_Table];             %#ok<AGROW>
    Summary_All     = [Summary_All; Summary_Table];         %#ok<AGROW>
    G_Shape_All     = [G_Shape_All; G_Shape_Table];         %#ok<AGROW>
    G_TimeScale_All = [G_TimeScale_All; G_TimeScale_Table]; %#ok<AGROW>
    Kappa_By_TTM{ttm_idx} = Kappa_Table;

    if isempty(R_Grid_All)
        R_Grid_All = R_Grid_Table;
    else
        if height(R_Grid_All) ~= height(R_Grid_Table)
            error('The R-grid length differs across TTM.');
        end
        grid_error = max(abs(R_Grid_All.gross_return_R - ...
            R_Grid_Table.gross_return_R));
        if grid_error > 1e-12
            error(['The 30,000-point R grid differs across TTM. ' ...
                'Maximum absolute difference = %.3e.'], grid_error);
        end
    end

    clear Data theta_hat kappa_vec SDF_Cell Model_Info
    clear Theta_Table Kappa_Table Summary_Table R_Grid_Table
    clear G_Shape_Table G_TimeScale_Table
end


%% Write only the compact handoff files

Kappa_All = combine_kappa_tables(Kappa_By_TTM);

writetable(Theta_All, fullfile(Path_Output, 'Theta_By_TTM.csv'));
writetable(Kappa_All, fullfile(Path_Output, 'Kappa_By_TTM.csv'));
writetable(Summary_All, fullfile(Path_Output, ...
    'CubicBSpline_EstimationSummary_AllTTM.csv'));
writetable(R_Grid_All, fullfile(Path_Output, ...
    'CubicBSpline_R_Grid.csv'));
writetable(G_Shape_All, fullfile(Path_Output, ...
    'CubicBSpline_GShape_By_TTM.csv'));
writetable(G_TimeScale_All, fullfile(Path_Output, ...
    'CubicBSpline_GTimeScale_By_TTM.csv'));

fprintf('\nAll four TTM estimations and compact CSV exports are complete.\n');
fprintf('Output folder: %s\n', Path_Output);
fprintf('No annual SDF-grid CSV files were generated.\n');


%% Local helper

function Wide_Table = combine_kappa_tables(Kappa_Tables)
%COMBINE_KAPPA_TABLES Keep the established wide format with unequal samples.

    max_rows = max(cellfun(@height, Kappa_Tables));
    Wide_Table = table();

    for table_idx = 1:numel(Kappa_Tables)
        Current_Table = Kappa_Tables{table_idx};
        variable_names = Current_Table.Properties.VariableNames;

        for variable_idx = 1:numel(variable_names)
            variable_name = variable_names{variable_idx};
            values = double(Current_Table.(variable_name));
            padded_values = NaN(max_rows, 1);
            padded_values(1:numel(values)) = values(:);
            Wide_Table.(variable_name) = padded_values;
        end
    end
end
