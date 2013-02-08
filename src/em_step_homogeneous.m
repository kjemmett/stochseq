function [S gamma xi LpX] = em_step_homogeneous(x, S0, A, e)
% [S g xi LpX] = em_step_homogeneous(x, S0, A, e)
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
%   gamma : marginal posterior matrix (T x L)
%   xi : joint posterior matrix (T x T x L)
%   LpX : log-likelihood (double)
%
% last modified : 2012-01-20


L = size(S0, 1); % dna sequence length
T = length(x); % read sequence length

% pre-calculate probabilities 
for d = 1:4
   M{d} = A * sparse(diag((1-e) * S0(:,d) + e/4));
end

% normalization constant
c = ones(T, 1); 

%%%%%%%%%%
% E Step %
%%%%%%%%%%

% Forward Pass
a = zeros(T, L);
a(1, 1) = 1;
for t = 2:T
    a(t, :) = a(t-1, :) * M{x(t)}; % recursion
    c(t) = sum(a(t, :)); % compute scaling factor
    a(t, :) = a(t, :) / c(t); % normalize alpha
    a(t,isnan(a(t,:))==1) = 0; % fix NaN
end

% Backward Pass
b = zeros(T, L); 
b(end, end) = 1; 
for t = T-1:-1:1
    b(t, :) = M{x(t+1)} * b(t+1,:)'; % recursion
    b(t, :) = b(t, :) / c(t+1); % normalize beta
    b(t,isnan(b(t,:))==1) = 0; % fix NaN
end

% gamma: marginal posterior
gamma = a .* b;
gamma = gamma / sum(gamma(1,:));
gamma(isnan(gamma)==1) = 0;

% xi: joint posterior
xi = zeros(L,L,T);
%for t=2:T
%    xitmp = c(t) * repmat(a(t-1,:)',1,L).*A.*repmat(S0(:,x(t))'.*b(t,:),L,1);
%    xi(:,:,t)=xitmp / sum(sum(xitmp));
%end
%xi(isnan(xi)==1) = 0;

% log likelihood
LpX = sum(log(c));

%%%%%%%%%%
% M Step %
%%%%%%%%%%

% update emission matrix
for d = 1:4
    S(:, d) = (sum(gamma(x==d, :),1) ./ sum(gamma, 1))';
end
S = S ./ repmat(sum(S,2),1,4);
