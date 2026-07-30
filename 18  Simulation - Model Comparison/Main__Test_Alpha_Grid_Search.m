clear; clc; close all;

%% Test purpose
% Compare a fixed 17-point alpha grid with a coarse-to-fine search,
% using the same data, DGP, B-spline specification, monotonicity
% restriction, and MLE routine as the formal simulation.


%% Paths

Path_MainFolder = ['D:\Google\我的雲端硬碟\學術｜研究與論文\論文著作\', ...
    'MLE Pricing Kernel'];

Path_Data = fullfile(Path_MainFolder, 'Code', '00  Output');
Path_RND = fullfile(Path_MainFolder, 'Code', '01  Output');
Path_Code_18 = fullfile(Path_MainFolder, 'Code', ...
    '18  Simulation - Model Comparison');
Path_Formal_Output = fullfile(Path_MainFolder, 'Code', '18  Output');
Path_Test_Output = fullfile(Path_MainFolder, 'Code', ...
    '18  Output - test');

addpath(Path_Code_18);

if ~exist(Path_Test_Output, 'dir')
    mkdir(Path_Test_Output);
end


%% Test settings

Target_TTM = 30;
b = 6;

Num_Test_Simulations = 3;
Test_Base_Seed = 42;

Fixed_Beta = 1.00;
Impose_Monotonicity = true;

Fixed_Alpha_Grid = round((0.60:0.05:1.40)', 10);
Coarse_Alpha_Grid = round((0.60:0.10:1.40)', 10);

assert(any(Fixed_Alpha_Grid == 1), ...
    'The fixed alpha grid must contain alpha = 1.');
assert(any(Coarse_Alpha_Grid == 1), ...
    'The coarse alpha grid must contain alpha = 1.');
assert(all(ismember(Coarse_Alpha_Grid, Fixed_Alpha_Grid)), ...
    'The coarse grid must be a subset of the fixed grid.');


%% Load realized returns

FileName = ['Realized_Return_TTM_', num2str(Target_TTM), '.csv'];
Realized_Return = readtable(fullfile(Path_Data, FileName));


%% Load risk-free gross factors

Risk_Free_Rate_All = readtable(fullfile( ...
    Path_Data, 'Risk_Free_GrossFactor_ByTargetTTM.csv'));


%% Load Q-measure PDF tables

Smooth_AllR = [];
Smooth_AllR_RND = [];

for year = 1996:2025
    input_filename = fullfile(Path_RND, sprintf( ...
        'TTM_%d_RND_Tables_%d.mat', Target_TTM, year));

    if exist(input_filename, 'file')
        Loaded = load(input_filename);
        Smooth_AllR = [Smooth_AllR, ...
            Loaded.Table_Smooth_AllR]; %#ok<AGROW>
        Smooth_AllR_RND = [Smooth_AllR_RND, ...
            Loaded.Table_Smooth_AllR_RND]; %#ok<AGROW>
    else
        warning('File %s does not exist.', input_filename);
    end
end

clear Loaded FileName input_filename year


%% Align RND, realized returns, and risk-free rates

rnd_var_names = string(Smooth_AllR_RND.Properties.VariableNames);
rnd_dates = str2double(regexp( ...
    rnd_var_names, '\d+', 'match', 'once'));

realized_dates = double(Realized_Return.date);

rf_date_col = ['date_', num2str(Target_TTM)];
rf_rate_col = ['rf_gross_TTM', num2str(Target_TTM)];

rf_dates_all = double(Risk_Free_Rate_All.(rf_date_col));
rf_values_all = double(Risk_Free_Rate_All.(rf_rate_col));

idx_rf_valid = isfinite(rf_dates_all) & isfinite(rf_values_all);
rf_dates_all = rf_dates_all(idx_rf_valid);
rf_values_all = rf_values_all(idx_rf_valid);

idx_master = ...
    ismember(realized_dates, rnd_dates) & ...
    ismember(realized_dates, rf_dates_all);

common_dates = sort(realized_dates(idx_master));

if isempty(common_dates)
    error('No common dates across the three input datasets.');
end

[tf_rnd, idx_rnd] = ismember(common_dates, rnd_dates);
[tf_ret, idx_ret] = ismember(common_dates, realized_dates);
[tf_rf, idx_rf] = ismember(common_dates, rf_dates_all);

if ~all(tf_rnd) || ~all(tf_ret) || ~all(tf_rf)
    error('Date alignment failed.');
end

Smooth_AllR = Smooth_AllR(:, idx_rnd);
Smooth_AllR_RND = Smooth_AllR_RND(:, idx_rnd);
Realized_Return = Realized_Return(idx_ret, :);
Risk_Free_Rate = rf_values_all(idx_rf);

R_observed = Realized_Return.realized_ret(:);
Risk_Free_Rate = Risk_Free_Rate(:);

fprintf('\nAligned sample: %d months, %08.0f to %08.0f.\n', ...
    numel(common_dates), common_dates(1), common_dates(end));

clear idx_master idx_rf_valid idx_rnd idx_ret idx_rf
clear tf_rnd tf_ret tf_rf rnd_dates realized_dates
clear rf_dates_all rf_values_all rf_date_col rf_rate_col
clear Risk_Free_Rate_All rnd_var_names


%% Global return range and precomputed B-spline basis

Global_Min_R = Inf;
Global_Max_R = -Inf;
months = Smooth_AllR.Properties.VariableNames;

for t = 1:numel(months)
    R_axis_t = Smooth_AllR.(months{t});
    Global_Min_R = min(Global_Min_R, min(R_axis_t));
    Global_Max_R = max(Global_Max_R, max(R_axis_t));
end

Global_Min_R = Global_Min_R * 0.9;
Global_Max_R = Global_Max_R * 1.1;

n_degree = 5;
k_order = n_degree + 1;
num_basis_function = b + 1;
num_breaks = num_basis_function - k_order + 2;
breaks = linspace(Global_Min_R, Global_Max_R, num_breaks);
knots = augknt(breaks, k_order);

T = numel(R_observed);
Basis_Precomputed = cell(T, 1);

for t = 1:T
    R_axis_t = Smooth_AllR.(months{t});
    Basis_Precomputed{t} = spcol( ...
        knots, k_order, R_axis_t(:));
end


%% Obtain the formal constrained-undistorted fitted null model

Observed_Fits_File = fullfile( ...
    Path_Formal_Output, 'Observed_Model_Fits.mat');

Loaded_Null_Fit = false;

if exist(Observed_Fits_File, 'file')
    Saved_Observed = load(Observed_Fits_File);

    if isfield(Saved_Observed, 'Observed_Fits') && ...
            isfield(Saved_Observed.Observed_Fits, 'Constraint') && ...
            isfield(Saved_Observed, 'Observed_Fit_Settings')

        Null_Fit = Saved_Observed.Observed_Fits.Constraint;
        Old_Settings = Saved_Observed.Observed_Fit_Settings;

        Loaded_Null_Fit = ...
            numel(Null_Fit.theta) == b + 1 && ...
            Null_Fit.alpha == 1 && ...
            Null_Fit.beta == 1 && ...
            Old_Settings.Target_TTM == Target_TTM && ...
            Old_Settings.b == b && ...
            isequal(Old_Settings.Sample_Dates(:), common_dates(:));
    end
end

if Loaded_Null_Fit
    theta_null = Null_Fit.theta;
    LL_null_observed = Null_Fit.LL;
    fprintf('Reusing the compatible formal constrained null fit.\n');
else
    fprintf('Compatible formal null fit not found; estimating it now.\n');

    [theta_null, LL_null_observed] = MLE_BSpline_estimation( ...
        Smooth_AllR, Smooth_AllR_RND, ...
        Realized_Return, Risk_Free_Rate, ...
        b, 1.00, 1.00, ...
        Global_Min_R, Global_Max_R, ...
        true, Basis_Precomputed);
end


%% Construct the common constrained-undistorted DGP

[R_Axis_Null, F_Null] = build_physical_cdf_bspline( ...
    theta_null, Risk_Free_Rate, Basis_Precomputed, ...
    Smooth_AllR, Smooth_AllR_RND, months, 1.00, 1.00);


%% Build the observed and simulated samples

Num_Datasets = 1 + Num_Test_Simulations;
Dataset_Names = strings(Num_Datasets, 1);
Return_Samples = cell(Num_Datasets, 1);

Dataset_Names(1) = "Observed";
Return_Samples{1} = R_observed;

for simulation_id = 1:Num_Test_Simulations
    rng(Test_Base_Seed + simulation_id, 'twister');

    Dataset_Names(simulation_id + 1) = sprintf( ...
        'Simulation_%03d', simulation_id);
    Return_Samples{simulation_id + 1} = ...
        draw_returns_from_monthly_cdf(R_Axis_Null, F_Null);
end


%% Inputs shared by the parallel workers

Search_Data = struct();
Search_Data.Smooth_AllR = Smooth_AllR;
Search_Data.Smooth_AllR_RND = Smooth_AllR_RND;
Search_Data.Realized_Return_Template = Realized_Return;
Search_Data.Risk_Free_Rate = Risk_Free_Rate;
Search_Data.Basis_Precomputed = Basis_Precomputed;
Search_Data.b = b;
Search_Data.Global_Min_R = Global_Min_R;
Search_Data.Global_Max_R = Global_Max_R;
Search_Data.Fixed_Beta = Fixed_Beta;
Search_Data.Impose_Monotonicity = Impose_Monotonicity;

if isempty(gcp('nocreate'))
    parpool;
end

Search_Data_Const = parallel.pool.Constant(Search_Data);
Progress_Queue = parallel.pool.DataQueue;

Completed_Datasets = 0;
Progress_Timer = tic;

afterEach(Progress_Queue, @show_dataset_progress);

fprintf(['\nStarting alpha-grid comparison for %d datasets ', ...
    '(observed + %d simulations).\n'], ...
    Num_Datasets, Num_Test_Simulations);
fprintf('Fixed grid: %d candidates.\n', numel(Fixed_Alpha_Grid));
fprintf('Coarse grid: %d candidates plus local midpoints.\n\n', ...
    numel(Coarse_Alpha_Grid));

Search_Results = cell(Num_Datasets, 1);

parfor dataset_idx = 1:Num_Datasets
    D = Search_Data_Const.Value;
    R_now = Return_Samples{dataset_idx};

    Simulated_Return = D.Realized_Return_Template;
    Simulated_Return.realized_ret = R_now;

    Fixed_Result = run_alpha_grid( ...
        Fixed_Alpha_Grid, D, Simulated_Return);

    % Reconstruct the exact coarse-to-fine decision from the fixed-grid
    % likelihoods. This avoids estimating identical alpha points twice.
    Adaptive_Result = derive_coarse_to_fine_result( ...
        Fixed_Result, Coarse_Alpha_Grid);

    Fixed_Seconds = sum(Fixed_Result.CandidateSeconds);
    Adaptive_Seconds = sum(Adaptive_Result.CandidateSeconds);

    One = struct();
    One.Dataset = Dataset_Names(dataset_idx);
    One.Fixed = Fixed_Result;
    One.Adaptive = Adaptive_Result;
    One.Fixed_Seconds = Fixed_Seconds;
    One.Adaptive_Seconds = Adaptive_Seconds;
    One.Alpha_Match = abs( ...
        Fixed_Result.SelectedAlpha - ...
        Adaptive_Result.SelectedAlpha) < 1e-10;
    One.Alpha_Difference = ...
        Adaptive_Result.SelectedAlpha - ...
        Fixed_Result.SelectedAlpha;
    One.LL_Difference = ...
        Adaptive_Result.SelectedLL - ...
        Fixed_Result.SelectedLL;

    Search_Results{dataset_idx} = One;
    send(Progress_Queue, dataset_idx);
end

delete(Search_Data_Const);


%% Build and export the comparison summary

Dataset = strings(Num_Datasets, 1);
Fixed_Selected_Alpha = nan(Num_Datasets, 1);
Adaptive_Selected_Alpha = nan(Num_Datasets, 1);
Alpha_Match = false(Num_Datasets, 1);
Alpha_Difference = nan(Num_Datasets, 1);
Fixed_Selected_LL = nan(Num_Datasets, 1);
Adaptive_Selected_LL = nan(Num_Datasets, 1);
Adaptive_Minus_Fixed_LL = nan(Num_Datasets, 1);
Fixed_Candidates = nan(Num_Datasets, 1);
Adaptive_Candidates = nan(Num_Datasets, 1);
Fixed_Minutes = nan(Num_Datasets, 1);
Adaptive_Estimated_Minutes = nan(Num_Datasets, 1);
Estimated_Speedup = nan(Num_Datasets, 1);
Adaptive_Selected_At_Search_Boundary = false(Num_Datasets, 1);

Details_Cells = cell(2 * Num_Datasets, 1);

for dataset_idx = 1:Num_Datasets
    One = Search_Results{dataset_idx};

    Dataset(dataset_idx) = One.Dataset;
    Fixed_Selected_Alpha(dataset_idx) = ...
        One.Fixed.SelectedAlpha;
    Adaptive_Selected_Alpha(dataset_idx) = ...
        One.Adaptive.SelectedAlpha;
    Alpha_Match(dataset_idx) = One.Alpha_Match;
    Alpha_Difference(dataset_idx) = One.Alpha_Difference;
    Fixed_Selected_LL(dataset_idx) = One.Fixed.SelectedLL;
    Adaptive_Selected_LL(dataset_idx) = One.Adaptive.SelectedLL;
    Adaptive_Minus_Fixed_LL(dataset_idx) = One.LL_Difference;
    Fixed_Candidates(dataset_idx) = numel(One.Fixed.AlphaGrid);
    Adaptive_Candidates(dataset_idx) = ...
        numel(One.Adaptive.AlphaGrid);
    Fixed_Minutes(dataset_idx) = One.Fixed_Seconds / 60;
    Adaptive_Estimated_Minutes(dataset_idx) = ...
        One.Adaptive_Seconds / 60;
    Estimated_Speedup(dataset_idx) = ...
        One.Fixed_Seconds / One.Adaptive_Seconds;
    Adaptive_Selected_At_Search_Boundary(dataset_idx) = ...
        One.Adaptive.SelectedAtSearchBoundary;

    Fixed_Details = table( ...
        repmat(One.Dataset, numel(One.Fixed.AlphaGrid), 1), ...
        repmat("Fixed_17_point", ...
            numel(One.Fixed.AlphaGrid), 1), ...
        One.Fixed.AlphaGrid, ...
        One.Fixed.LogLikelihoodGrid, ...
        One.Fixed.ExitFlagGrid, ...
        One.Fixed.CandidateSeconds, ...
        'VariableNames', { ...
            'Dataset', 'Method', 'Alpha', ...
            'LogLikelihood', 'ExitFlag', 'CandidateSeconds'});

    Adaptive_Details = table( ...
        repmat(One.Dataset, numel(One.Adaptive.AlphaGrid), 1), ...
        repmat("Coarse_to_fine", ...
            numel(One.Adaptive.AlphaGrid), 1), ...
        One.Adaptive.AlphaGrid, ...
        One.Adaptive.LogLikelihoodGrid, ...
        One.Adaptive.ExitFlagGrid, ...
        One.Adaptive.CandidateSeconds, ...
        'VariableNames', { ...
            'Dataset', 'Method', 'Alpha', ...
            'LogLikelihood', 'ExitFlag', 'CandidateSeconds'});

    Details_Cells{2 * dataset_idx - 1} = Fixed_Details;
    Details_Cells{2 * dataset_idx} = Adaptive_Details;
end

Candidate_Details = vertcat(Details_Cells{:});

Summary = table( ...
    Dataset, ...
    Fixed_Selected_Alpha, Adaptive_Selected_Alpha, ...
    Alpha_Match, Alpha_Difference, ...
    Fixed_Selected_LL, Adaptive_Selected_LL, ...
    Adaptive_Minus_Fixed_LL, ...
    Fixed_Candidates, Adaptive_Candidates, ...
    Fixed_Minutes, Adaptive_Estimated_Minutes, ...
    Estimated_Speedup, ...
    Adaptive_Selected_At_Search_Boundary);

Summary_File = fullfile(Path_Test_Output, ...
    'Alpha_Search_Comparison_Summary.csv');
Details_File = fullfile(Path_Test_Output, ...
    'Alpha_Search_Candidate_Details.csv');
MAT_File = fullfile(Path_Test_Output, ...
    'Alpha_Search_Test_Results.mat');

writetable(Summary, Summary_File);
writetable(Candidate_Details, Details_File);

Test_Settings = struct( ...
    'Target_TTM', Target_TTM, ...
    'b', b, ...
    'Num_Test_Simulations', Num_Test_Simulations, ...
    'Test_Base_Seed', Test_Base_Seed, ...
    'Fixed_Beta', Fixed_Beta, ...
    'Impose_Monotonicity', Impose_Monotonicity, ...
    'Fixed_Alpha_Grid', Fixed_Alpha_Grid, ...
    'Coarse_Alpha_Grid', Coarse_Alpha_Grid, ...
    'Sample_Dates', common_dates(:), ...
    'Null_LogLikelihood_Observed', LL_null_observed);

save(MAT_File, ...
    'Summary', 'Candidate_Details', 'Search_Results', ...
    'Test_Settings', '-v7.3');

fprintf('\nAlpha-grid comparison completed.\n\n');
disp(Summary);
fprintf('Summary CSV:\n%s\n\n', Summary_File);
fprintf('Candidate details CSV:\n%s\n\n', Details_File);
fprintf('MAT results:\n%s\n', MAT_File);


%% Local functions

function Result = derive_coarse_to_fine_result( ...
    Fixed_Result, coarse_grid)

    [is_coarse, coarse_idx] = ismember( ...
        coarse_grid, Fixed_Result.AlphaGrid);

    if ~all(is_coarse)
        error('Every coarse-grid point must exist in the fixed grid.');
    end

    coarse_ll = Fixed_Result.LogLikelihoodGrid(coarse_idx);
    coarse_exit = Fixed_Result.ExitFlagGrid(coarse_idx);
    coarse_valid = isfinite(coarse_ll) & coarse_exit > 0;

    if ~any(coarse_valid)
        error('No valid coarse-grid estimates were obtained.');
    end

    coarse_selection_ll = coarse_ll;
    coarse_selection_ll(~coarse_valid) = -Inf;
    [~, coarse_best_idx] = max(coarse_selection_ll);
    coarse_best = coarse_grid(coarse_best_idx);
    coarse_step = median(diff(coarse_grid));
    fine_step = coarse_step / 2;

    refinement_grid = [ ...
        coarse_best - fine_step; ...
        coarse_best + fine_step];

    refinement_grid = round(refinement_grid, 10);
    refinement_grid = refinement_grid( ...
        refinement_grid >= min(coarse_grid) & ...
        refinement_grid <= max(coarse_grid));
    refinement_grid = setdiff( ...
        refinement_grid, coarse_grid, 'stable');

    adaptive_alpha = sort(unique([coarse_grid; refinement_grid]));
    [is_adaptive, adaptive_idx] = ismember( ...
        adaptive_alpha, Fixed_Result.AlphaGrid);

    if ~all(is_adaptive)
        error('Adaptive-grid point missing from the fixed grid.');
    end

    all_alpha = Fixed_Result.AlphaGrid(adaptive_idx);
    all_ll = Fixed_Result.LogLikelihoodGrid(adaptive_idx);
    all_exit = Fixed_Result.ExitFlagGrid(adaptive_idx);
    all_seconds = Fixed_Result.CandidateSeconds(adaptive_idx);

    valid = isfinite(all_ll) & all_exit > 0;

    if ~any(valid)
        selected_alpha = NaN;
        selected_ll = -Inf;
    else
        selection_ll = all_ll;
        selection_ll(~valid) = -Inf;
        [selected_ll, best_idx] = max(selection_ll);
        selected_alpha = all_alpha(best_idx);
    end

    Result = struct();
    Result.AlphaGrid = all_alpha;
    Result.LogLikelihoodGrid = all_ll;
    Result.ExitFlagGrid = all_exit;
    Result.CandidateSeconds = all_seconds;
    Result.SelectedAlpha = selected_alpha;
    Result.SelectedLL = selected_ll;
    Result.CoarseSelectedAlpha = coarse_best;
    Result.SelectedAtSearchBoundary = ...
        selected_alpha == min(coarse_grid) || ...
        selected_alpha == max(coarse_grid);
end


function Result = run_alpha_grid(alpha_grid, D, Return_Table)

    alpha_grid = alpha_grid(:);
    num_alpha = numel(alpha_grid);

    ll_grid = -Inf(num_alpha, 1);
    exit_grid = -999 * ones(num_alpha, 1);
    candidate_seconds = nan(num_alpha, 1);

    for alpha_idx = 1:num_alpha
        alpha_now = alpha_grid(alpha_idx);

        Candidate_Timer = tic;

        [~, ll_now, ~, exit_now] = ...
            MLE_BSpline_estimation( ...
            D.Smooth_AllR, D.Smooth_AllR_RND, ...
            Return_Table, D.Risk_Free_Rate, ...
            D.b, alpha_now, D.Fixed_Beta, ...
            D.Global_Min_R, D.Global_Max_R, ...
            D.Impose_Monotonicity, ...
            D.Basis_Precomputed);

        candidate_seconds(alpha_idx) = toc(Candidate_Timer);

        ll_grid(alpha_idx) = ll_now;
        exit_grid(alpha_idx) = exit_now;
    end

    valid = isfinite(ll_grid) & exit_grid > 0;

    if ~any(valid)
        selected_alpha = NaN;
        selected_ll = -Inf;
    else
        selection_ll = ll_grid;
        selection_ll(~valid) = -Inf;
        [selected_ll, best_idx] = max(selection_ll);
        selected_alpha = alpha_grid(best_idx);
    end

    Result = struct();
    Result.AlphaGrid = alpha_grid;
    Result.LogLikelihoodGrid = ll_grid;
    Result.ExitFlagGrid = exit_grid;
    Result.CandidateSeconds = candidate_seconds;
    Result.SelectedAlpha = selected_alpha;
    Result.SelectedLL = selected_ll;
end


function R_sim = draw_returns_from_monthly_cdf( ...
    R_axis_cell, F_cell)

    T = numel(R_axis_cell);
    R_sim = zeros(T, 1);

    for t = 1:T
        R_axis_t = R_axis_cell{t}(:);
        F_t = F_cell{t}(:);

        [F_unique, idx_unique] = unique(F_t, 'stable');
        R_unique = R_axis_t(idx_unique);

        if numel(F_unique) < 2
            error([ ...
                'Physical CDF at t=%d has fewer than ', ...
                'two unique values.'], t);
        end

        u = rand;
        R_sim(t) = interp1( ...
            F_unique, R_unique, u, 'linear', 'extrap');
    end
end


function show_dataset_progress(~)

    persistent Completed_Count
    persistent Progress_Timer

    if isempty(Completed_Count)
        Completed_Count = 0;
        Progress_Timer = tic;
    end

    Completed_Count = Completed_Count + 1;

    fprintf(['Alpha-search dataset completed: %d; ', ...
        'elapsed = %.2f minutes.\n'], ...
        Completed_Count, toc(Progress_Timer) / 60);
end
