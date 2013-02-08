function S = h_step(em)
% S = h_step(em)
%
% average sequence inference of each read
%
% inputs: 
%   em : struct
%
% outputs:
%   S : averaged emission matrix (Lx4)
%
% last modified : 2012-01-20

S = mean(cat(3, em(:).S), 3);
