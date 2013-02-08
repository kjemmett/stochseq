NEWEST = zeros(100,4);

for i=1:5

    NEWEST = NEWEST + build_seqest(max_gamma(inf5.em(i).gamma),model.reads(i).x);

end

NEWEST = NEWEST ./ repmat(sum(NEWEST,2),1,4);
