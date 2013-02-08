clear all


L=10;
p=0.70;
e=0.0;
N=2;

dna = [2 3 2 3 2 1 2 4 2 1];
model = stochseq_build(L,p,e,N,'dna',dna);

A = full(gen_transmatrix(L,p));

S = 0.25*ones(4,L);
S(:,1) = [0 0 0 0];
S(:,2) = [0 0 0 0];
S(:,end) = [0 0 0 0];
S(:,end-1) = [0 0 0 0];
S(model.dna(1),1) = 1;
S(model.dna(2),2) = 1;
S(model.dna(end),end) = 1;
S(model.dna(end-1),end-1) = 1;

% store length of each read sequence (should really use an absorbing state here
% to normalize the lengths)
for i=1:N
    T(i) = length(model.reads(i).z);
    a{i} = zeros(T(i), L);
    a{i}(1,1) = 1;
    a{i}(2,2) = 1;
    c{i} = ones(T(i),1);
    b{i} = zeros(T(i), L);
    b{i}(end,end) = 1;
end

% precompute probabilities
for d=1:4
    M{d} = A * diag((1-e) * S(d,:) + e/4);
end

% forward pass for each read
for i=1:N
    for t=3:T(i)
        a{i}(t,:) = a{i}(t-1,:) * M{model.reads(i).x(t)};
        c{i}(t) = sum(a{i}(t,:));
        a{i}(t,:) = a{i}(t,:) / c{i}(t);
    end
end

for i=1:N
    for t = T(i)-1:-1:1
        b{i}(t, :) = M{model.reads(i).x(t+1)} * b{i}(t+1,:)';
        b{i}(t, :) = b{i}(t,:) / c{i}(t+1);
    end
end

for i=1:N
    gamma{i} = a{i} .* b{i} ;
    gamma{i} = gamma{i} / sum(gamma{i}(1,:));
end

for i=1:N
    xi{i} = zeros(L,L,T(i));
    ds{i} = zeros(T(i),2);
    for t=2:T(i)

        xitmp = c{i}(t) * repmat(a{i}(t-1,:)',1,L) .* A .* ... 
                repmat(S(model.reads(i).x(t),:) .* b{i}(t,:),L,1);

        xi{i}(:,:,t) = xitmp / sum(sum(xitmp));

        ds{i}(t,:) = sum([diag(xi{i}(:,:,t),1) diag(xi{i}(:,:,t),-1)]);
    end
end
