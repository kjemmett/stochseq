proj_path = '/Users/kje/columbia/research/wigginslab/stochseq/src';
addpath(genpath(proj_path))

N_range = [1 2 5 10 25 50 75 100];

p = 0.75;
e = 0.05;
L = 100;

% generate dna for parameter sweep
dna = gen_dna(L);

% preallocate memory for matrices
for i=1:length(N_range)
    for j=1:5
        model = stochseq_build(L,p,e,N_range(i),'dna',dna);
        inf = stochseq_infer(model,'epsilon',1e-3);
        sweep(i,j).N = N_range(i);
        sweep(i,j).ed = inf.Edit_Distance;
        clear model,inf
    end
end

