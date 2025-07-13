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
Risk_Free_Rate = readtable(fullfile(Path_Data_01_main, FileName));

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

clear FileName input_filename year Path_Data_01 Path_Data_01_main Path_Data_02 data


%% Interpolate f^*_t at each realized gross return R_{t+1}

T = height(Realized_Return);
R_realized_vec = Realized_Return.realized_ret;     % Realized gross returns R_{t+1}
date_vec = Realized_Return.date;                   % Date vector

f_star_interp = zeros(T, 1);                       % Interpolated f^*_t(R_{t+1})

for t = 1:T
    date_str = num2str(date_vec(t));
    R_grid = Smooth_AllR.(date_str);               % Grid of R (R_{t+1} variable, not realized)
    f_star_grid = Smooth_AllR_RND.(date_str);      % Corresponding f^*_t(R) values

    % Ensure R_realized is within interpolation range
    R_realized = R_realized_vec(t);
    if R_realized < min(R_grid) || R_realized > max(R_grid)
        warning('t = %d: R_{t+1} = %.4f is out of bounds. Assigned tiny density.', t, R_realized);
        f_star_interp(t) = 1e-10;
    else
        f_star_interp(t) = interp1(R_grid, f_star_grid, R_realized, 'linear');
    end
end

clear t T R_realized R_grid f_star_grid date_str


%% Estimate θ = (c_1, ..., c_N) by maximizing log-likelihood

% Add paths
Path_Code_02 = fullfile(Path_MainFolder, 'Code', '02  Maximum Likelihood Estimation');
addpath(Path_Code_02);

% Prepare input for estimation
R_vec   = Realized_Return.realized_ret;                                    % R_{t+1}
Rf_vec  = Risk_Free_Rate.rate;                                             % R_f^t
dates   = Realized_Return.date;                                            % date t
N       = 3;                                                               % Degree of polynomial (can adjust)

% Initial guess for θ: [c_1, ..., c_N]
theta0 = zeros(N, 1);

% Set optimization options
options = optimoptions('fmincon', ...
    'Display', 'iter-detailed', ...
    'Algorithm', 'interior-point', ...
    'SpecifyObjectiveGradient', false);

% Define bounds if needed (optional)
LB = [];    % e.g., -Inf(N+1, 1)
UB = [];    % e.g., +Inf(N+1, 1)

% Objective function (note: fmincon minimizes, so we use negative LL)
obj_fun = @(theta) -log_likelihood_function(theta, R_vec, Rf_vec, ...
                        f_star_interp, N, Smooth_AllR, Smooth_AllR_RND, dates);

% Run fmincon
[theta_hat, neg_LL, exitflag, output] = fmincon(obj_fun, theta0, [], [], [], [], LB, UB, [], options);

% Report result
disp('Estimated parameters (θ):');
disp(theta_hat);
disp(['Final log-likelihood = ', num2str(-neg_LL)]);

% Optional: save results
% save('MLE_Result.mat', 'theta_hat', 'neg_LL');
