function inf = playiterate2(model)
% function inf_iter = playiterate2(model)

L = model.seqlength;
p = model.bias;
e = model.err;
N = model.nreads;
dna = model.dna;

reads = model.reads;


S = 0.25 * ones(L,4);
S(1,:) = [0 0 0 0];
S(end,:) = [0 0 0 0];
S(1,dna(1)) = 1;
S(end,dna(end)) = 1;

A = gen_transmatrix(L,p);

for i=1:N
    T(i) = length(reads(i).z);
    alpha{i} = zeros(T(i),L);
    alpha{i}(1,1) = 1;
    alpha{i}(2,2) = 2;
end

entropy(1) = L + 100;
epsilon = .001;
q=1;
for t=2:max(T)
    for j=2:10
        % compute alpha
        for i=1:N
            if t<=T(i)
                alpha{i}(t,:) = alpha{i}(t-1,:) * A * diag(S(:,reads(i).x(t)));
                c{i}(t) = sum(alpha{i}(t,:));
                alpha{i}(t,:) = alpha{i}(t,:) / c{i}(t);
            end
        end

        % update S
        for i=1:N
            if t<=T(i)
                S(:,reads(i).x(t)) = S(:,reads(i).x(t)) + alpha{i}(t,:)';
            end
        end
        %S = S ./ repmat(sum(S,2),1,4);
        inf.h(q).S = S;
        inf.h(q).inf_ent = calc_entropy(S);
        q=q+1;
        % update alpha
        for i=1:N
            if t<=T(i)
                alpha{i}(t,:) = alpha{i}(t-1,:) * A * diag(S(:,reads(i).x(t)));
                c{i}(t) = sum(alpha{i}(t,:));
                alpha{i}(t,:) = alpha{i}(t,:) / c{i}(t);
            end
        end
        % convergence criteria : sequence entropy
        entropy(j) = sum(calc_entropy(S));
        fprintf('at t= %d, j= %d, entropy= %e\n',t,j,entropy(j));
        %if abs(entropy(j)-entropy(j-1)) < epsilon * abs(entropy(j-1))
        %    break
        %end
    end
end

inf.T = T;
inf.S = S;
inf.alpha = alpha;
