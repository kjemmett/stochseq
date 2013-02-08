function [z] = build_pathest(tv)
%
% build pathest vector from transition vector

T = length(tv)+1;
L = sum(tv)+1;

z = zeros(T,1);
z(1) = 1;

zm = zeros(T,L);
zm(1,1) = 1;

for i=2:T
    z(i) = z(i-1) + tv(i-1);
    zm(i,z(i)) = 1;
end
