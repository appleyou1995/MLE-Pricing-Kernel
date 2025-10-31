%% Main Function: MLE gamma Estimation for Pricing Kernel

function [gamma_hat, log_lik, delta_vec, M_vec] = MLE_gamma_estimation( ...
    Smooth_AllR, Smooth_AllR_RND, Realized_Return, Risk_Free_Rate, ...
    L, use_delta, alpha, beta)

    % Settings
    rng(0);
    dates  = Realized_Return.date;
    R_vec  = Realized_Return.realized_ret;
    Rf_vec = Risk_Free_Rate;

    % Initial guess
    gamma0 = zeros(L,1);
    LB = [];
    UB = [];

    % Optimization options
    options = optimoptions('fmincon', ...
        'Display', 'iter-detailed', ...
        'Algorithm', 'interior-point', ...
        'SpecifyObjectiveGradient', false);

    % Define objective function (fmincon minimizes, so we negate the log-likelihood)
    obj_fun = @(param) -log_likelihood_function(param, R_vec, Rf_vec, L, ...
        Smooth_AllR, Smooth_AllR_RND, dates, use_delta, ...
        alpha, beta);

    % Run optimization
    [gamma_hat, neg_LL, exitflag, ~] = fmincon(obj_fun, gamma0, [], [], [], [], LB, UB, [], options);

    % Return log-likelihood (positive value)
    log_lik = -neg_LL;

    % Check optimization result
    if     exitflag > 0,  disp('✅ fmincon converged successfully.');
    elseif exitflag == 0, disp('⚠️ fmincon reached the iteration limit.');
    else,                 disp(['❌ fmincon failed. Exit flag: ', num2str(exitflag)]);
    end
    
    % Display output
    disp('Estimated coefficients:');
    disp(gamma_hat);
    disp(['Final log-likelihood = ', num2str(log_lik)]);

    % Post-estimation call
    [~, delta_vec, M_vec] = log_likelihood_function(gamma_hat, R_vec, Rf_vec, L, ...
    Smooth_AllR, Smooth_AllR_RND, dates, use_delta, alpha, beta);

end