function S = build_seqest(zm,x)
% function S = build_seqest(zm,x)
%
% build a sequence estimate from an inferred path matrix

[T L] = size(zm);

S = zeros(L,4);

for d = 1:4
    S(:,d) = (sum(zm(x==d,:),1) ./ sum(zm,1))';
end

%S = S ./ repmat(sum(S,2),1,4);
