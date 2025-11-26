function [risk_pref, deriv] = compute_risk_pref_from_M(R_axis, M_bar)

    % Make sure they are column vectors
    R_axis = R_axis(:);
    M_bar  = M_bar(:);

    % ---- 1. Numerical derivatives using finite differences ----
    M1 = gradient(M_bar, R_axis);   % dM/dR
    M2 = gradient(M1,    R_axis);   % d^2 M / dR^2
    M3 = gradient(M2,    R_axis);   % d^3 M / dR^3

    % ---- 2. Relative thresholds to avoid division by near-zero derivatives ----
    % Threshold for M (rarely triggered)
    epsM  = 1e-12;
    M_bar_safe = M_bar;
    % M_bar_safe(abs(M_bar) < epsM) = NaN;

    % Threshold for M' based on its scale
    M1_safe = M1;
    % M1_safe(abs(M1) < epsM) = NaN;

    % Threshold for M''
    M2_safe = M2;
    % M2_safe(abs(M2) < epsM) = NaN;

    % ---- 3. Compute risk preference measures ----
    % ARA(R) = - M'(R) / M(R)
    ARA = - M1_safe ./ M_bar_safe;

    % RRA(R) = - R * M'(R) / M(R)
    RRA = - R_axis .* (M1_safe ./ M_bar_safe);

    % AP(R) = - M''(R) / M'(R)
    AP  = - M2_safe ./ M1_safe;

    % RP(R) = - R * M''(R) / M'(R)
    RP  = - R_axis .* (M2_safe ./ M1_safe);

    % AT(R) = - M'''(R) / M''(R)
    AT  = - M3 ./ M2_safe;

    % RT(R) = - R * M'''(R) / M''(R)
    RT  = - R_axis .* (M3 ./ M2_safe);

    % ---- 4. Remove unstable endpoints ----
    % idx_end = [1; length(R_axis)];
    % ARA(idx_end) = NaN;
    % RRA(idx_end) = NaN;
    % AP(idx_end)  = NaN;
    % RP(idx_end)  = NaN;
    % AT(idx_end)  = NaN;
    % RT(idx_end)  = NaN;

    % ---- 5. Output ----
    risk_pref = struct();
    risk_pref.R_axis = R_axis;
    risk_pref.ARA = ARA;
    risk_pref.RRA = RRA;
    risk_pref.AP  = AP;
    risk_pref.RP  = RP;
    risk_pref.AT  = AT;
    risk_pref.RT  = RT;

    deriv = struct();
    deriv.M1 = M1;
    deriv.M2 = M2;
    deriv.M3 = M3;
end
