proj_path = '/Users/kje/columbia/research/wigginslab/stochseq'
addpath(genpath(proj_path))

L=1000;

p=0.8;
e=0.05;

N = 5;

model = stochseq_build(L,p,e,N);
inf = stochseq_infer(model);

