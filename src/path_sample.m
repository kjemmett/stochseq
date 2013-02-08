function path_sample(gamma,x,dna,)

[T L] = size(gamma);

% compute most likely path from gamma
[ignore path] = max(gamma,[],2);

% construct transition vector from path
for i=2:length(path)
    trans(i) = path(i)-path(i-1);
end

% construct path matrix from most likely path
new_gamma = zeros(T,L);
for i=1:T
    new_gamma(i,path(i)) = 1;

% now take x and compute sequence estimate
seq = zeros(L,4);
for i=1:length(path)
    seq(path(i),x(i)) = seq(path(i),x(i)) + 1;
end

seq = seq ./ repmat(sum(seq,2),1,4);

