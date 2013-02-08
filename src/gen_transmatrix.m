function A = gen_transmatrix(L, p)
% A = gen_transmatrix(L, p)
%
% Generates transition matrix
%
% inputs:
%   L : sequence length
%   p : forward bias
%
% outputs:
%   A : transition matrix (LxL)
%
% last modified : 2012-01-20

A = sparse(toeplitz([0 1-p zeros(1, L-2)], [0 p zeros(1, L-2)]));

% move forward with prob 1 from start
A(1, 1:2) = [0 1]; 

% move backward with prob 1-p from end
A(L, L-1:L) = [1-p 0]; 
