function [c, ceq] = monotone_nonlcon(gamma, L)

    % 1) Define the checking interval: x = log(R) ∈ [log(0.8), log(1.2)]
    a = log(0.8);
    b = log(1.2);

    % 2) Choose the number of grid points K:
    %    Larger K gives a denser check but adds more constraints and computation
    K = 20;
    x = linspace(a, b, K).';   % K×1 vector, each element is a grid point x_k

    % 3) Compute the derivative of the polynomial:
    %    poly_sum'(x) = Σ_{l=1}^L l * gamma_l * x^(l-1)
    dpoly = zeros(K,1);
    for l = 1:L
        dpoly = dpoly + l * gamma(l) .* (x.^(l-1));
    end

    % 4) Set a small tolerance to avoid infeasibility due to numerical precision
    tol = 1e-6;

    % 5) fmincon requires c(x) ≤ 0.
    %    We want dpoly - tol ≥ 0, which is equivalent to -(dpoly - tol) ≤ 0.
    c = -(dpoly - tol);   % K×1 vector of inequality constraints
    ceq = [];             % No equality constraints

end
