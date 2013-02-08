function dna_matrix = gen_dnamatrix(dna)
% dna_matrix = gen_dnamatrix(dna)
%
% Takes an input dna sequence and constructs a
% binary output matrix representing the input
% sequence.
%
% inputs:
%	dna : (Lx1) int
%		DNA sequence. List of integers from 1
%		to 4 representing input base pairs
%
% outputs:
%	dna_matrix : (Lx4) binary
%		DNA sequence represented as binary
%		matrix. Four columns representing four
%		base pairs, one nonzero element per row.
%

L = size(dna,1);
dna_matrix = zeros(L,4);

for d=1:4
    dna_matrix(dna(:)==d,d)=1;
end
