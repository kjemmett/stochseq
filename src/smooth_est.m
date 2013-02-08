function S_new = smooth_est(S,e)
% function S_new = smooth_est(S,e)
%
% input:
%   S : Lx4 sequence estimate
%   e : smoothing factor

L = size(S,1);

noise = (1-e)/4 * ones(size(S));

S_temp = S*e + noise;

S_new = zeros(size(S));

S_new(1,:) = (S_temp(1,:) + S_temp(2,:))/2;

for i = 2:L-1
    S_new(i,:) = (S_temp(i-1,:) + S_temp(i,:) + S_temp(i+1,:))/3;
end

S_new(L,:) = (S_temp(L-1,:) + S_temp(L,:))/2;
