function DebugTable = make_debug_table(all_results, L_target, alpha_target, beta_target, x_min, x_max)
% MAKE_DEBUG_TABLE
%   Extracts derivatives and risk preference measures for a specific
%   (L, alpha, beta) combination over a given R-interval [x_min, x_max].
%
%   INPUT:
%       all_results : cell array produced previously (each cell contains:
%                     L, alpha, beta, R_axis, M_bar, risk_pref, deriv)
%       L_target    : target L (integer)
%       alpha_target: target alpha (double)
%       beta_target : target beta (double)
%       x_min       : lower bound of R (double)
%       x_max       : upper bound of R (double)
%
%   OUTPUT:
%       DebugTable : table containing R, M, M1, M2, M3, ARA, RRA, AP, RP, AT, RT
%
%   USAGE:
%       T = make_debug_table(all_results, 3, 1.00, 0.90, 1.15, 1.20);

    % ---- 1. Find the matching entry ----
    idx_found = [];
    for idx = 1:numel(all_results)
        res = all_results{idx};
        if res.L == L_target && res.alpha == alpha_target && res.beta == beta_target
            idx_found = idx;
            break;
        end
    end

    if isempty(idx_found)
        error('Case not found: L=%d alpha=%.2f beta=%.2f', ...
            L_target, alpha_target, beta_target);
    end

    res = all_results{idx_found};

    % ---- 2. Extract interval mask ----
    R = res.R_axis(:);
    mask = (R >= x_min) & (R <= x_max);

    % ---- 3. Extract values ----
    R_sub  = R(mask);
    M_sub  = res.M_bar(:); M_sub = M_sub(mask);

    M1_sub = res.deriv.M1(:); M1_sub = M1_sub(mask);
    M2_sub = res.deriv.M2(:); M2_sub = M2_sub(mask);
    M3_sub = res.deriv.M3(:); M3_sub = M3_sub(mask);

    rp = res.risk_pref;

    ARA_sub = rp.ARA(:); ARA_sub = ARA_sub(mask);
    RRA_sub = rp.RRA(:); RRA_sub = RRA_sub(mask);
    AP_sub  = rp.AP(:);  AP_sub  = AP_sub(mask);
    RP_sub  = rp.RP(:);  RP_sub  = RP_sub(mask);
    AT_sub  = rp.AT(:);  AT_sub  = AT_sub(mask);
    RT_sub  = rp.RT(:);  RT_sub  = RT_sub(mask);

    % ---- 4. Assemble table ----
    DebugTable = table( ...
        R_sub, M_sub, ...
        M1_sub, M2_sub, M3_sub, ...
        ARA_sub, RRA_sub, AP_sub, RP_sub, AT_sub, RT_sub, ...
        'VariableNames', { ...
            'R', 'M', ...
            'M1', 'M2', 'M3', ...
            'ARA', 'RRA', 'AP', 'RP', 'AT', 'RT'});

end
