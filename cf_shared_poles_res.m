function [zi, ci_list, r_inf] = cf_shared_poles_res(fct_list, n, varargin)
% CF_SHARED_POLES_RES  Shared-denominator CF on (-inf,0] with z-plane LS residues.
%   [zi, ci_list, r_inf] = cf_shared_poles_res(fct_list, n, ...
%                       'K', 75, 'nf', 1024, 'sigma', 9, 'weights', [], 'fitN', 1e4)
%
% Inputs
%   fct_list : cell {f0, f1, ..., fp}. Each f accepts vector x<=0.
%   n        : denominator degree (number of shared poles).
%
%   'K'      : # Chebyshev coeffs used for CF pole extraction 
%   'nf'     : FFT grid size on unit circle for CF pole extraction 
%   'sigma'  : scale in x = sigma*(t-1)/(t+1)       
%   'weights': (p+1)-vector to weight stacked Hankels 
%   'fitN'   : # of real-axis points for z-plane LS 
%
% Outputs
%   zi       : [n x 1] shared poles in the z-plane 
%   ci_list  : n-by-(p+1) matrix; ci_list(:,j) is [n x 1] residues for f_j
%   r_inf    : (p+1)-by-1 constants (value at +inf) for each f_j
%
%   Reference:
%
%   Awad H. Al-Mohy, Shared-Pole Carath\'eodory--Fej\'er Approximations 
%           for Linear Combinations of φ-Functions, 2025.
%
%=======================================================================

% ---------------- options -----------------
K       = getn(varargin,   'K'      ,   100);
nf      = getn(varargin,   'nf'     ,  2^10);
sigma   = getn(varargin, 'sigma'    ,     9);
fitN    = getn(varargin,   'fitN'   ,   1e4);
m       = numel(fct_list);
wts     = getn(varargin, 'weights'  , ones(m,1));
wts     = wts(:);

% -------- sanity checks --------
if nf <= K
    error('cf_shared_realaxis:badGrid', ...
        'Require nf > K (got nf=%d, K=%d). Increase nf or decrease K.', nf, K);
end
if n >= K
    error('cf_shared_realaxis:badDegree', ...
        'Degree n must be < K (got n=%d, K=%d). Increase K or reduce n.', n, K);
end
if numel(wts) ~= m || any(~isfinite(wts)) || any(wts<=0)
    error('cf_shared_realaxis:badWeights', ...
        'weights must be a positive vector of length numel(fct_list).');
end

% -------- (A) CF pole extraction: map, sample, stacked Hankel, SVD -----
% Unit-circle samples and Chebyshev abscissae (twice over)
w = exp(2i*pi*(0:nf-1)/nf);     % |w| = 1, nf roots of unity
g = (real(w) ~= -1);
x = zeros(size(w));
x(g) = sigma*real(( (w(g)-1)./(w(g)+1) ).^2);     % |w| = 1 --> (-inf,0]

% Chebyshev coefficients via FFT
C = zeros(K+1, m);
for j = 1:m
    F = zeros(size(w));
    Fj = feval(fct_list{j}, x(g), varargin{:});    
    F(g) = Fj;
    c = real(fft(F))/nf;            % c_0..c_{nf-1}
    C(:,j) = c(1:K+1);              % keep first K+1
end

% Build stacked Hankel 
Hstack = zeros(m*K, K);
for j = 1:m
    cj = C(:,j);
    Hj = hankel(cj(2:K+1));          % KxK Hankel from c_1..c_K
    rows = (j-1)*K + (1:K);
    Hstack(rows, :) = sqrt(wts(j)) * Hj;
end

% SVD on tall-skinny Hstack
[~,~,V] = svd(Hstack, 'econ');   

v = V(:, n+1).';                    % polynomial coeffs (row), desc. powers

% Poles on the circle: zeros with |q|>1; choose n nearest to unit circle
q_all = roots(v);                   % length K-1
q_out = q_all(abs(q_all) > 1);
if numel(q_out) < n
    error('cf_shared_pole_res:fewOutsideRoots', ...
        ['Only %d zeros with |q|>1 found; need n=%d....' ...
        ' Increase K/nf or adjust sigma.'], numel(q_out), n);
end
[~,ord] = sort(abs(q_out), 'ascend');  % closeness to unit circle
qj = q_out(ord(1:n));

% Map to z-plane poles (right half-plane) and order by imag 
zi = sigma*( (qj-1)./(qj+1) ).^2;       
[~, order] = sort(imag(zi), 'descend');
zi = zi(order);

% -------- (B) Residues via z-plane least squares on (-inf,0] ----------
ci_list = zeros(n,m);
r_inf   = zeros(m,1);

% Real-axis fitting grid: dense near 0 and -inf (symmetric log-spacing)
half = max(2, floor(fitN/2));
C = 100;                        % user choice, e.g. 50 or 100
M = C * sigma;                  % far-left magnitude proportional to sigma
m1 = sigma * logspace(-12, 0, half);     % near-zero band, scaled
m2 = sigma * logspace(0, log10(M/sigma), half);  % 1 .. M/sigma
x_fit = [-m1(end:-1:1), -m2].';
x_fit = unique([x_fit; 0], 'stable');  % drop duplicate -1 and append 0
L = numel(x_fit);

% Design matrix A = [1, 1/(x-zi_1), ..., 1/(x-zi_n)]
W = 1 ./ (x_fit - zi(:).');            % L x n
A = [ones(L,1), W];                    % L x (n+1)
[Q, R] = qr(A,0);                      % thin QR for repeated use
for j = 1:m
    fj = feval(fct_list{j}, x_fit, varargin{:});    
   % coeff = A\ fj;                    % complex LS
    coeff         = R\(Q'*fj); 
    r_inf(j)      = coeff(1);
    ci_list(:,j)  = coeff(2:end);
end

% Enforce conjugacy on each pair by averaging
Ipos  = find(imag(zi) >  0);
Ineg  = find(imag(zi) <  0);

% Build a map from positive-imag poles to their conjugate partners
neg_pool = Ineg(:).';
for jp = Ipos(:).'
    [~,idx] = min(abs(zi(jp) - conj(zi(neg_pool))));
    mate(jp) = neg_pool(idx);              %#ok<AGROW>
    neg_pool(idx) = [];                    % consume once
end

for j = 1:m
    eta = ci_list(:,j);    
    for jp = Ipos(:).'
        jm = mate(jp);
        eta_avg = 0.5*(eta(jp) + conj(eta(jm)));
        eta(jp) = eta_avg;
        eta(jm) = conj(eta_avg);
    end
    ci_list(:,j) = eta;
end

end

% --------------- helper ----------------
function val = getn(vargs, name, default)
idx = find(strcmpi(vargs(1:2:end), name), 1);
if isempty(idx), val = default; else, val = vargs{2*idx}; end
end
