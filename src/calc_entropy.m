function entropy = calc_entropy(est)
%function entropy = calc_entropy(est)
%
% Calculate the normalized entropy of an input sequence estimate
%
% Inputs
% ------
%
% est : (L x 4)
%	Sequence estimate of length L.
%
% Output
% ------
%
% entropy: (1 x L)
%	Normalized entropy
%
% last modified : 2012-01-20

L = size(est,1); % input est should be Lx4
entropy = zeros(1,L);
q=est.*log(est);
q(isnan(q)==1)=0; % necessary to fix (0,1,0,0) elements
entropy = -sum(q,2)/log(4); % divide by log(4) to normalize result
