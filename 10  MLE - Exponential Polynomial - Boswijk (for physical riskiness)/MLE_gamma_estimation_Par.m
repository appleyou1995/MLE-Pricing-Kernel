%% Main Function: MLE gamma Estimation for Pricing Kernel
function [gamma_hat, log_lik, kappa_vec, M_cell] = MLE_gamma_estimation_Par( ...
    Smooth_AllR, Smooth_AllR_RND, Realized_Return, Risk_Free_Rate, L)

    % Settings
    rng(0);
    dates  = Realized_Return.date;
    R_vec  = Realized_Return.realized_ret;
    Rf_vec = Risk_Free_Rate;
    T      = length(R_vec);

    % --- Data Pre-processing for parfor ---
    disp('Converting Struct to Cell Array for parallel computing...');
    R_grids_All = cell(T, 1);
    f_star_All  = cell(T, 1);
    
    for t = 1:T
        date_str = num2str(dates(t));
        R_grids_All{t} = Smooth_AllR.(date_str);
        f_star_All{t}  = Smooth_AllR_RND.(date_str);
    end
    
    % Initial guess
    gamma0 = zeros(L,1);
    LB = [];
    UB = [];

    % Optimization options
    options = optimoptions('fmincon', ...
        'Display', 'iter-detailed', ...
        'Algorithm', 'interior-point', ...
        'SpecifyObjectiveGradient', false);

    % Define objective function
    obj_fun = @(param) -log_likelihood_function_par(param, R_vec, Rf_vec, L, ...
        R_grids_All, f_star_All);

    % Run optimization
    disp('Starting fmincon with parallel computing...');
    [gamma_hat, neg_LL, exitflag, ~] = fmincon(obj_fun, gamma0, [], [], [], [], LB, UB, [], options);

    % Return log-likelihood
    log_lik = -neg_LL;

    % Check optimization result
    if     exitflag > 0,  disp('fmincon converged successfully.');
    elseif exitflag == 0, disp('fmincon reached the iteration limit.');
    else,                 disp(['fmincon failed. Exit flag: ', num2str(exitflag)]);
    end
    
    % Display output
    disp('Estimated coefficients:');
    disp(gamma_hat);
    disp(['Final log-likelihood = ', num2str(log_lik)]);

    % Post-estimation call
    [~, kappa_vec, M_cell] = log_likelihood_function_par(gamma_hat, R_vec, Rf_vec, L, ...
        R_grids_All, f_star_All);
end