proj_path = '/Users/kje/columbia/research/wigginslab/stochseq'
addpath(genpath(proj_path))

L = 10;
e = 0.0;

p_range = 0.55:0.05:1.0;
N_range = [1 2 4 6 8 10];

dna = gen_dna(L);

for i=1:length(p_range)
    for j=1:length(N_range)
        model = stochseq_build(L,p_range(i),e,N_range(j),'dna',dna);
        for k=1:10
            fprintf('p = %d, N = %d, k = %d\n',p_range(i),N_range(j),k);
            inf = stochseq_infer(model,'epsilon',1e-3,'debug',0);
            ed_array(i,j,k) = inf.Edit_Distance;
            fprintf('edit distance = %d\n',ed_array(i,j,k));
        end
    end
end
