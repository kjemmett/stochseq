function [S, g, LpX] = em_step_triplet(x, S0, A, e)
% [S g LpX] = em_step_triplet(x, S0, A, e)
%
% executes a single EM step on a single DNA read sequence.
% uses a stationary transition matrix
%
% inputs:
%   x  : observed state sequence (T x 1)
%   S0 : initial emission matrix (L x 4)
%   A  : transition matrix (L x L)
%   e  : base-call error rate (double)
%
% outputs:
%   S : updated emission matrix (L x 4)
%   g : marginal posterior matrix (T x L)
%   xi : joint posterior matrix (T x T x L)
%   LpX : log-likelihood (double)
%
% last modified : 2013-03-04

L = size(S0, 1); % dna sequence length
T = length(x);   % read sequence length

% initialize matrix with observation probabilities
px_z = S0(:, x(:, 2))';
%pxp_z = S0(:, x(1:end-1, 3))';
%pxm_z = S0(:, x(2:end, 1))';
%px_z(1:end-1, :) = (px_z(1:end-1, :) + S0(:, x(2:end, 1))') / 2;
%px_z(2:end, :) = px_z(2:end, :) + S0(:, x(1:end-1, 3))' / 3;

%px_z = ((1-e) * S0(:, x(:,2)) + e/4)';
%px_z(1:end-1, :) = px_z(1:end-1, :) + ((1-e) * S0(:, x(2:end, 1)) + e/4)';
%px_z(2:end, :) = px_z(2:end, :) + ((1-e) * S0(:, x(1:end-1, 3)) + e/4)';


% set boundary conditions
a0 = zeros(L, 1);
a0(1) = 1;
bT = zeros(L, 1);
bT(end) = 1;

% run forward backward
%[g, xi, LpX] = forwback_banded(px_z, A, [-1 1], a0, bT);
[g, LpX] = forwback_sparse(px_z, A, A', a0, bT);
%[g, LpX] = forwaback_stochseq(px_z, p);

% update emission matrix
for d = 1:4
    S(:, d) = (sum(g(x(:, 2)==d, :), 1) ./ sum(g, 1))';

    tmp = (sum(g(x(:, 3)==d, :), 1) ./ sum(g, 1))';
    S(2:end, d) = S(2:end, d) + tmp(1:end-1);
    tmp = (sum(g(x(:, 1)==d, :), 1) ./ sum(g, 1))';
    S(1:end-1, d) = S(1:end-1, d) + tmp(2:end);
end
S = S ./ repmat(sum(S,2),1,4);
