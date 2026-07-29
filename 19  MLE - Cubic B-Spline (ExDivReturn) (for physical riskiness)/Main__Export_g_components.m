clear; clc;

%% Paths and settings

Path_MainFolder = ['D:\Google\我的雲端硬碟\學術｜研究與論文\論文著作\' ...
    'MLE Pricing Kernel'];
Path_Data = fullfile(Path_MainFolder, 'Code', '00  Output');
Path_SDF_Output = fullfile(Path_MainFolder, 'Code', '19  Output');
Path_Output = fullfile(Path_MainFolder, 'Code', '19  Output - g');

TTM_List = [30, 60, 90, 180];
Expected_Grid_Size = 30000;
Spline_Degree = 3;

if ~exist(Path_Output, 'dir'), mkdir(Path_Output); end


%% Shared input and output containers

Risk_Free_All = readtable(fullfile( ...
    Path_Data, 'Risk_Free_GrossFactor_ByTargetTTM.csv'));

R_Grid_Table = table();
G_Shape_All = table();
G_TimeScale_All = table();

fprintf('Exporting compact g components for TTM = %s\n', ...
    mat2str(TTM_List));


%% Construct the return-only shape and monthly time scale

for ttm_idx = 1:numel(TTM_List)
    Target_TTM = TTM_List(ttm_idx);
    prefix = sprintf('CubicBSpline_TTM%03d', Target_TTM);

    fprintf('\nTTM = %d\n', Target_TTM);

    Theta_Table = readtable(fullfile( ...
        Path_SDF_Output, [prefix '_Theta.csv']));
    Knot_Table = readtable(fullfile( ...
        Path_SDF_Output, [prefix '_Knots.csv']));
    Kappa_Table = readtable(fullfile( ...
        Path_SDF_Output, [prefix '_Kappa.csv']));
    Manifest_Table = readtable(fullfile( ...
        Path_SDF_Output, [prefix '_SDFGrid_Manifest.csv']));

    require_columns(Theta_Table, {'TTM', 'theta'}, ...
        [prefix '_Theta.csv']);
    require_columns(Knot_Table, {'knot_index', 'knot_value'}, ...
        [prefix '_Knots.csv']);
    require_columns(Manifest_Table, {'calendar_year', 'relative_file'}, ...
        [prefix '_SDFGrid_Manifest.csv']);

    theta = double(Theta_Table.theta(:));
    Knot_Table = sortrows(Knot_Table, 'knot_index');
    knots = double(Knot_Table.knot_value(:)');
    spline_order = numel(knots) - numel(theta);

    if spline_order ~= Spline_Degree + 1
        error(['TTM = %d implies spline order %d; a cubic B-spline ' ...
            'requires order 4.'], Target_TTM, spline_order);
    end

    % Read one existing annual SDF file to recover and verify the common
    % 30,000-point gross-return grid used by the estimation output.
    Manifest_Table = sortrows(Manifest_Table, 'calendar_year');
    relative_grid_file = char(string(Manifest_Table.relative_file(1)));
    Sample_Grid = readtable(fullfile(Path_SDF_Output, relative_grid_file));
    require_columns(Sample_Grid, ...
        {'quote_date', 'grid_index', 'gross_return_R', 'sdf_M'}, ...
        relative_grid_file);

    sample_quote_date = double(Sample_Grid.quote_date(1));
    sample_rows = double(Sample_Grid.quote_date) == sample_quote_date;
    sample_grid_index = double(Sample_Grid.grid_index(sample_rows));
    sample_R = double(Sample_Grid.gross_return_R(sample_rows));
    sample_M = double(Sample_Grid.sdf_M(sample_rows));

    [sample_grid_index, sort_order] = sort(sample_grid_index);
    sample_R = sample_R(sort_order);
    sample_M = sample_M(sort_order);

    if numel(sample_R) ~= Expected_Grid_Size || ...
            ~isequal(sample_grid_index, (1:Expected_Grid_Size)')
        error('Unexpected R grid for TTM = %d.', Target_TTM);
    end
    if any(~isfinite(sample_R)) || any(diff(sample_R) <= 0)
        error('The R grid must be finite and strictly increasing.');
    end

    if isempty(R_Grid_Table)
        R_Grid_Table = table(sample_grid_index, sample_R, ...
            'VariableNames', {'grid_index', 'gross_return_R'});
    else
        grid_error = max(abs(R_Grid_Table.gross_return_R - sample_R));
        if grid_error > 1e-12
            error('The R grid differs across TTM; maximum error = %.3e.', ...
                grid_error);
        end
    end

    % q(R) = sum_i theta_i B_i(R). Use the B-form spline and its
    % analytical derivatives to avoid numerical-difference noise.
    spline_curve = spmak(knots, theta');
    q = fnval(spline_curve, sample_R')';
    q_d1 = fnval(fnder(spline_curve, 1), sample_R')';
    q_d2 = fnval(fnder(spline_curve, 2), sample_R')';

    % Confirm that spmak/fnval reproduces the basis representation used in
    % the MLE code: spcol(knots, order, R) * theta.
    q_basis = spcol(knots, spline_order, sample_R) * theta;
    basis_error = max(abs(q - q_basis));
    if basis_error > 1e-11
        error('B-spline representation mismatch for TTM = %d: %.3e.', ...
            Target_TTM, basis_error);
    end

    % h(R) = exp(-q(R)); g_t(R) = g_scale_t * h(R).
    g_shape = exp(-q);
    g_shape_d1 = -q_d1 .* g_shape;
    g_shape_d2 = (q_d1 .^ 2 - q_d2) .* g_shape;

    if any(~isfinite(g_shape)) || any(g_shape <= 0) || ...
            any(~isfinite(g_shape_d1)) || any(~isfinite(g_shape_d2))
        error('Non-finite g-shape values found for TTM = %d.', Target_TTM);
    end

    monotonicity_tolerance = 1e-10 * max(1, max(abs(g_shape_d1)));
    if min(g_shape_d1) < -monotonicity_tolerance
        error('g-shape is not non-decreasing for TTM = %d.', Target_TTM);
    end

    Shape_Table = table( ...
        repmat(Target_TTM, Expected_Grid_Size, 1), ...
        sample_grid_index, g_shape, g_shape_d1, g_shape_d2, ...
        'VariableNames', {'TTM', 'grid_index', 'g_shape', ...
        'g_shape_d1', 'g_shape_d2'});
    G_Shape_All = [G_Shape_All; Shape_Table]; %#ok<AGROW>

    date_col = sprintf('date_TTM_%d', Target_TTM);
    exdate_col = sprintf('exdate_TTM_%d', Target_TTM);
    kappa_col = sprintf('kappa_TTM_%d', Target_TTM);
    rf_date_col = sprintf('date_%d', Target_TTM);
    rf_col = sprintf('rf_gross_TTM%d', Target_TTM);

    require_columns(Kappa_Table, {date_col, exdate_col, kappa_col}, ...
        [prefix '_Kappa.csv']);
    require_columns(Risk_Free_All, {rf_date_col, rf_col}, ...
        'Risk_Free_GrossFactor_ByTargetTTM.csv');

    quote_date = double(Kappa_Table.(date_col));
    expiration_date = double(Kappa_Table.(exdate_col));
    kappa = double(Kappa_Table.(kappa_col));
    rf_date = double(Risk_Free_All.(rf_date_col));
    rf_value = double(Risk_Free_All.(rf_col));

    valid_rf = isfinite(rf_date) & isfinite(rf_value);
    rf_date = rf_date(valid_rf);
    rf_value = rf_value(valid_rf);

    if numel(unique(quote_date)) ~= numel(quote_date)
        error('Duplicate kappa quote dates found for TTM = %d.', Target_TTM);
    end
    if numel(unique(rf_date)) ~= numel(rf_date)
        error('Duplicate risk-free quote dates found for TTM = %d.', ...
            Target_TTM);
    end

    [date_is_present, rf_location] = ismember(quote_date, rf_date);
    if ~all(date_is_present)
        missing_dates = quote_date(~date_is_present);
        error('Missing risk-free factors for TTM = %d; first date = %08d.', ...
            Target_TTM, missing_dates(1));
    end

    rf_gross = rf_value(rf_location);
    if any(~isfinite(kappa)) || any(~isfinite(rf_gross)) || any(rf_gross <= 0)
        error('Invalid kappa or risk-free values for TTM = %d.', Target_TTM);
    end

    g_scale = exp(-kappa) ./ rf_gross;
    if any(~isfinite(g_scale)) || any(g_scale <= 0)
        error('Invalid g time scale for TTM = %d.', Target_TTM);
    end

    TimeScale_Table = table( ...
        repmat(Target_TTM, numel(quote_date), 1), ...
        quote_date, expiration_date, rf_gross, kappa, g_scale, ...
        'VariableNames', {'TTM', 'quote_date', 'expiration_date', ...
        'rf_gross', 'kappa', 'g_scale'});
    G_TimeScale_All = [G_TimeScale_All; TimeScale_Table]; %#ok<AGROW>

    % Verify the split representation against both its direct formula and
    % one existing monthly SDF grid.
    max_split_error = 0;
    for month_idx = 1:numel(quote_date)
        g_direct = 1 ./ (rf_gross(month_idx) .* ...
            exp(kappa(month_idx) + q));
        g_from_components = g_scale(month_idx) .* g_shape;
        max_split_error = max(max_split_error, ...
            max(abs(g_direct - g_from_components)));
    end

    sample_kappa_row = find(quote_date == sample_quote_date, 1);
    if isempty(sample_kappa_row)
        error('Sample SDF quote date is absent from kappa for TTM = %d.', ...
            Target_TTM);
    end
    reconstructed_M = exp(kappa(sample_kappa_row) + q);
    sdf_relative_error = max(abs(reconstructed_M - sample_M) ./ ...
        max(1, abs(sample_M)));

    if max_split_error > 1e-12 || sdf_relative_error > 1e-11
        error(['Validation failed for TTM = %d: split error %.3e, ' ...
            'SDF relative error %.3e.'], Target_TTM, ...
            max_split_error, sdf_relative_error);
    end

    fprintf(['  months = %d, grid points = %d, basis error = %.3e, ' ...
        'split error = %.3e, SDF relative error = %.3e\n'], ...
        numel(quote_date), numel(sample_R), basis_error, ...
        max_split_error, sdf_relative_error);
end


%% Export compact CSV files

writetable(R_Grid_Table, fullfile(Path_Output, ...
    'CubicBSpline_R_Grid.csv'));
writetable(G_Shape_All, fullfile(Path_Output, ...
    'CubicBSpline_GShape_By_TTM.csv'));
writetable(G_TimeScale_All, fullfile(Path_Output, ...
    'CubicBSpline_GTimeScale_By_TTM.csv'));

fprintf('\nSaved compact g components to:\n%s\n', Path_Output);
fprintf('  R grid rows: %d\n', height(R_Grid_Table));
fprintf('  g-shape rows: %d\n', height(G_Shape_All));
fprintf('  g-time-scale rows: %d\n', height(G_TimeScale_All));


%% Local helper

function require_columns(Input_Table, Required_Columns, File_Label)
    missing_columns = setdiff(Required_Columns, ...
        Input_Table.Properties.VariableNames);
    if ~isempty(missing_columns)
        error('%s is missing columns: %s', File_Label, ...
            strjoin(missing_columns, ', '));
    end
end
