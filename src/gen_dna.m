function dna = gen_dna(L)
% dna = gen_dna(L)
%
% Generates a random 'dna sequence', where each base is represented
% by an integer between 1 and 4.
%
% Inputs
% ------
%
% 	L : int
%		Sequence length.
%
% Outputs
% -------
%
%	dna : (L x 1) int
%   	Sequence. Each element is int from 1 to 4.

dna = ceil(4*rand(L, 1));