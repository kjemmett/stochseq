%proj_path = '/Users/kje/columbia/research/wigginslab/stochseq'
proj_path = '/proj/kje2109/stochseq'

addpath(genpath(proj_path))

fid = fopen([proj_path,'/scripts/nreads_sweep.txt'],'w');

L = 100;
e = 0.05;

prange = 0.75:0.05:1.0;
Nrange = [1 2 5 10 15 20 25 30];



dna = gen_dna(L);



for i=1:length(prange)
    for j=1:length(Nrange)
        for k=1:5
            fprintf('p=%e\t%N=%d\tk=%d\n',prange(i),Nrange(j),k);
            model = stochseq_build(L,prange(i),e,Nrange(j),'dna',dna);
            inf = stochseq_infer(model,'debug',0);
            fprintf(fid,'%e %n %e\n',prange(i),Nrange(j),inf.Edit_Distance/L);
        end
    end
end

fclose(fid);
