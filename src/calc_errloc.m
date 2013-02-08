function el = calc_errloc(S,dna)
% function el = calc_errloc(S,dna)
%
% Identifies locations of inference errors
%
% Inputs
% ------
%
% S : (Lx4)
%	Inferred DNA sequence
%
% dna : Lx1
%	True DNA sequence
%
% Outputs
% -------
%
% el : Lx1
% 	Locations of inference errors
%
% ec : int
%	
[ignore,guess]=max(S,[],2);
el = guess - dna;
el = el~=0;
sum(el)