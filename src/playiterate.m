reads = model.reads;
dna = model.dna;

A = full(gen_transmatrix(L,p));

for i=1:N
    T(i) = length(reads(i).z);
    alpha{i} = zeros(T(i),L);
    alpha{i}(1,1) = 1;
    alpha{i}(2,2) = 1;
    alpha2{i} = alpha{i};
    c{i} = ones(T(i),1);
end

S = 0.25*ones(4,L);
S(:,1) = 0 * S(:,1);
%S(:,1:2) = 0 * S(:,1:2);
S(:,end-1:end) = 0 * S(:,end-1:end);
S(dna(1),1) = 1;
%S(dna(2),2) = 1;
S(dna(end),end) = 1;
S(dna(end-1),end-1) = 1;

S2 = S;
S2(:,3:end-2) = zeros(4,L-4);

for t=2:max(T)
    % compute the alpha term
    for i=1:N
        if t<=T(i)
            alpha{i}(t,:) = alpha{i}(t-1,:) * A * diag(S(reads(i).x(t),:));
            c{i}(t) = sum(alpha{i}(t,:));
            alpha{i}(t,:) = alpha{i}(t,:) / c{i}(t);
        end
    end

    % use this to update S
    for i=1:N
        if t<=T(i)
            S(reads(i).x(t),:) = S(reads(i).x(t),:) + alpha{i}(t,:);
            S2(reads(i).x(t),:) = S2(reads(i).x(t),:) + alpha{i}(t,:);
        end
    end
    S = S ./ repmat(sum(S,1),4,1);
    S2 = S2 ./ repmat(sum(S2,1),4,1);
    S2(isnan(S2)) = 0;
    
    %plot_sequence_inference_static(S',dna)
    %pause

    % update alpha term
    for i=1:N
        if t<=T(i)
            alpha2{i}(t,:) = alpha{i}(t-1,:) * A * diag(S2(reads(i).x(t),:));
            c{i}(t) = sum(alpha2{i}(t,:));
            alpha2{i}(t,:) = alpha2{i}(t,:) / c{i}(t);
        end
    end

end

% need heuristic to prevent updating S term if not clear which direction traveled.

% backwards step


