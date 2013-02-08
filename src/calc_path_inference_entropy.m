function [path_inf_entropy] = calculate_path_inference_entropy(gamma)
% function [path_inf_entropy] = calculate_path_inference_entropy(gamma)
%
% Calculates normalized entropy in path inference from gamma matrix
%
% INPUTS
% ------
%
% gamma : (T x L) double
%
% Outputs
% -------
%
% path_inf_entropy : (T x 1) double
%
T = size(gamma,1);
L = size(gamma,2);
ent=gamma.*log(gamma);
ent(isnan(ent)==1) = 0;
path_inf_entropy = -sum(ent,2)/log(L);
