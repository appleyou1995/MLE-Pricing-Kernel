%% Compute 10-bin PIT statistics

function Stats = compute_pit_statistics(pit_vec, num_bins)

    if nargin < 2 || isempty(num_bins)
        num_bins = 10;
    end

    pit_vec = pit_vec(:);
    pit_vec = pit_vec(isfinite(pit_vec));

    if isempty(pit_vec)
        error('No finite PIT observations are available.');
    end

    pit_vec = min(max(pit_vec, 0), 1);

    edges = linspace(0, 1, num_bins + 1);
    counts = histcounts(pit_vec, edges);
    bin_probability = counts ./ numel(pit_vec);

    uniform_probability = 1 / num_bins;
    sup_uniform = max(abs( ...
        bin_probability - uniform_probability));

    % One-sample Kolmogorov-Smirnov statistic against Uniform(0,1).
    u = sort(pit_vec);
    n = numel(u);
    upper_gap = (1:n)' ./ n - u;
    lower_gap = u - (0:n-1)' ./ n;
    ks_stat = max([upper_gap; lower_gap]);

    Stats = struct();
    Stats.BinProbability = bin_probability;
    Stats.SupUniform     = sup_uniform;
    Stats.KS             = ks_stat;
    Stats.NumObservations = n;
end
