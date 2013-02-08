function [M,S] = make_pathsetmatrix(model)

N = size(model.reads,2);

lmax = 0;
for i=1:N
    l = length(model.reads(i).z);
    if l>lmax
        lmax = l;
    end
end

M = zeros(N,lmax);
S = zeros(N,lmax);

for i=1:N
    M(i,1:length(model.reads(i).z)) = model.reads(i).z';
    S(i,1:length(model.reads(i).x)) = model.reads(i).x';
end


