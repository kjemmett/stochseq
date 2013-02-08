function M = max_gamma(gamma)
% function M = max_gamma(gamma)
%
% makes a path_est matrix from
% an input gamma matrix

M = zeros(size(gamma));

for i = 1:size(M,1)
    [a b] = max(gamma(i,:),[],2);
    M(i,b) = 1;
end
