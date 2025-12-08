%CFPHIMV Carath\'eodory--Fej\'er Rational Approximations of φ-Functions.
%   Compute  w = sum_{j=0}^p φ_j(A) v_j
%   using shared-pole partial fractions.
%
% Inputs
%   phidx   : vector of the required indices of φ_j.
%   A       : n-by-n matrix (sparse/dense) with eigenvalues in the
%             left-half plane.
%   V       : n-by-(p+1)  matrix  [v0 v1 ... vp]  (each v_j is n-by-1)
%   nPoles  : number of poles (must be even for real data and <= 14)
%
% Output
%   w       : the n-by-1 vector w = sum_{j=0}^p φ_j(A) v_j.
%
% Example:
%      N  = 500; h = 1e2;
%      A  = -h*gallery('poisson', N);     % 2D Laplacian 
%      phidx = [0, 1, 3, 7];
%      V  = randn( size(A,1),numel(phidx) );
%
%      w = cfphimv(phidx, A, V);
%
%      O = zeros(size(A,1),1);
%      f = phimv(1,1,A,V(:,1),[V(:,2), O, V(:,3), O, O, O, V(:,end)]);
%      norm(w - f,1)/norm(f,1)
%
%      Both w and f approximate
%      φ_0(A)*V(:,1) + φ_1(A)*V(:,2) + φ_3(A)*V(:,3) + φ_7(A)*V(:,4)
%
% NOTE: phimv is available in https://github.com/aalmohy/phimv
%
% See also: PHI, CF_SHARED_POLES_RES functions.
%
%   Reference:
%
%   Awad H. Al-Mohy, Shared-Pole Carath\'eodory--Fej\'er Approximations 
%           for Linear Combinations of φ-Functions, 2025.
%
