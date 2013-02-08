stochseq: random walks on DNA
=============================

To run:

* ``stochseq_build(seqlength,bias,err,nreads)`` : builds
* ``stochseq_infer(model)`` : performs inference on model object

Functions:

* ``gen_dna(L)`` : generates a random dna sequence where each base is 
            represented by an integer between 1 and 4
* ``gen_read(dna, bias, err)`` : generates a random walk along 'dna'

Parameters:

* ``seqlength`` : sequence length
* ``bias`` : forward bias
* ``err`` : per base call error rate
* ``nreads`` : number of reads

TODO:

* Scale to L > 100 (goal is L=50000)
* Parallelize multiple read EM
* 3-mer window states
