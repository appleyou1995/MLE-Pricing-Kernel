%% Stable trapezoidal integration when the integrand is stored in logs

function log_integral = log_trapz_exp(x, log_y)

    x = x(:);
    log_y = log_y(:);

    if numel(x) ~= numel(log_y) || numel(x) < 2
        error('x and log_y must have the same length and at least two points.');
    end

    if any(isnan(log_y)) || any(log_y == Inf)
        error('log_y may contain -Inf for zero density, but not NaN or +Inf.');
    end

    dx = diff(x);

    if any(~isfinite(dx)) || any(dx <= 0)
        error('x must be finite and strictly increasing.');
    end

    left  = log_y(1:end-1);
    right = log_y(2:end);

    pair_max = max(left, right);
    log_pair_sum = -Inf(size(pair_max));

    finite_pair = isfinite(pair_max);

    log_pair_sum(finite_pair) = pair_max(finite_pair) + log( ...
        exp(left(finite_pair)  - pair_max(finite_pair)) + ...
        exp(right(finite_pair) - pair_max(finite_pair)));

    log_terms = log(dx / 2) + log_pair_sum;
    term_max = max(log_terms);

    if ~isfinite(term_max)
        log_integral = -Inf;
        return
    end

    log_integral = term_max + log(sum(exp(log_terms - term_max)));
end