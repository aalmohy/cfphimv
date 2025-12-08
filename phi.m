function y = phi(x, k)
%PHI  Robust evaluation of the scalar phi_k(x) function.
%   y = phi(x,k) for array x (real or complex) and integer k>=0.
%
%   phi_k(x) = sum_{j=0}^∞ x^j/(j+k)!.
%   Strategy:
%     - x == 0          : phi_k(0) = 1/k!
%     - |x| <= 0.5      : adaptive Taylor series
%     - |x| >  0.5      : expm1-based phi_1 + upward recurrence

    if ~(isscalar(k) && k >= 0 && k == floor(k))
        error('k must be a nonnegative integer scalar.');
    end

    % Ensure floating input (avoid integer-typed outputs).
    if ~isfloat(x), x = double(x); end

    y = zeros(size(x), 'like', x);

    % ---- exact zeros
    mask0 = (x == 0);
    if any(mask0(:))
        y(mask0) = cast(invfact(k), 'like', x);  % 1/k!
    end

    % ---- remaining points
    mask = ~mask0;
    if ~any(mask(:)), return; end

    % Build full-size masks for small/large based on |x|
    small_full = false(size(x));
    small_full(mask) = abs(x(mask)) <= 0.5;
    large_full = mask & ~small_full;

    % Small magnitude: adaptive Taylor
    if any(small_full(:))
        y(small_full) = phi_series_adaptive(x(small_full), k);
    end

    % Moderate/large magnitude: expm1 + upward recurrence
    if any(large_full(:))
        y(large_full) = phi_hybrid_stable(x(large_full), k);
    end
end

% ==================== helpers (local functions) ====================

function ys = phi_series_adaptive(x, k)
% Adaptive series: phi_k(x) = sum_{j>=0} x^j/(j+k)! with term recurrence.
    tol  = eps(cast(1,'like',x));                  % GPU/precision-safe
    ys   = cast(invfact(k), 'like', x) .* ones(size(x), 'like', x);  % j=0 term
    term = ys;

    j = 0;
    maxJ = 2000;  % safety cap
    active = true(size(x));
    while any(active(:)) && j < maxJ
        j = j + 1;
        % term_j = term_{j-1} .* (x / (k+j))
        term(active) = term(active) .* ( x(active) ./ (k + j) );     
        ys(active)   = ys(active) + term(active);
        active(active) = abs(term(active)) > tol .* (abs(ys(active)) + 1);
    end
    if j == maxJ
        warning('phi:seriesConvergence', ...
            'Taylor series for k=%d reached maximum iterations (%d).', k, maxJ);
    end
end

function yl = phi_hybrid_stable(x, k)
% For |x| > 0.5: phi_1 via expm1/x (with tiny-x series), then upward recurrence.
    if k == 0
        yl = exp(x);
        return
    end

    % --- phi_1 (stable)
    % Use expm1(x)./x when |x| is not tiny; otherwise use short series.
    phi1 = ones(size(x), 'like', x);
    nz = abs(x) > 1e-8;
    if any(nz(:))
        phi1(nz) = expm1(x(nz)) ./ x(nz);
    end
    if any(~nz(:))
        xs = x(~nz);
        % 1 + x/2 + x^2/6 + x^3/24 (enough for ~1e-24 when |x|<=1e-8)
        phi1(~nz) = 1 + xs./2 + (xs.^2)/6 + (xs.^3)/24;
    end
    if k == 1
        yl = phi1;
        return
    end

    % --- Upward recurrence:
    %     phi_j = (phi_{j-1} - 1/(j-1)!)/x,  j = 2..k
    invF = invfact_vec(k, x);  % invF(j+1) = 1/j!, typed like x
    prev = phi1;
    % (x == 0) cannot occur here because this branch is only for |x|>0.5
    for j = 2:k
        prev = (prev - invF(j)) ./ x;   % invF(j) = 1/(j-1)!
    end
    yl = prev;
end

function v = invfact_vec(k, proto)
% Return [1/0!, 1/1!, ..., 1/k!] safely, cast like proto.
    v = cast(exp(-gammaln((0:k) + 1)), 'like', proto);
end

function a = invfact(q)
% Return 1/q! safely (double); caller may cast.
    a = exp(-gammaln(q + 1));
end
