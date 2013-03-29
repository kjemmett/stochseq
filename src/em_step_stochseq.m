function [S, g, LpX] = em_step_sparse(x, S0, p, e)
% [S g LpX] = em_step_sparse(x, S0, p, e)
%
% executes a single EM step on a single DNA read sequence.
% uses a stationary transition matrix
%
% inputs:
%   x  : observed state sequence (T x 1)
%   S0 : initial emission matrix (L x 4)
%   p  : forward transition probability (1 x 1)
%   e  : base-call error rate (double)
%
% outputs:
%   S : updated emission matrix (L x 4)
%   g : marginal posterior matrix (T x L)
%   xi : joint posterior matrix (T x T x L)
%   LpX : log-likelihood (double)
%
% last modified : 2012-01-20

L = size(S0, 1); % dna sequence length
T = length(x); % read sequence length

% initialize matrix with observation probabilities
px_z = ((1-e) * S0(:, x) + e/4)';

% run forward backward
[g, LpX] = forwback_stochseq(px_z, p);

% update emission matrix
for d = 1:4
    S(:, d) = (sum(g(x==d, :),1) ./ sum(g, 1))';
end
S = S ./ repmat(sum(S,2),1,4);
