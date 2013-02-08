function [model] = stochseq_build(seqlength,bias,err,nreads,varargin)
% function [model] = stochseq_build(seqlength,bias,err,nreads,varargin)
%
% data generation for stochseq
%  1. generates random dna sequence (or uses input sequence)
%  2. builds a set of 'nreads' dna sequences
%
% use output struct as input to stochseq_infer
%
% inputs:
%   seqlength : sequence length
%   bias : forward bias
%   err : error rate
%   nreads : number of reads to generate
%
% varargin:
%   'dna' : skip gen_dna step, use input dna sequence
%   'debug' : binary, controls stdout messages [default = 1]
%
% outputs:
%   model : struct containing read sequences and data parameters
%       model.reads : struct containing read sequence information
%           model.reads.x : observed state sequence
%           model.reads.z : latent state sequence
%           model.reads.e : read error locations
%           model.reads.t : read sequence transition vector
%       model.dna : dna sequence vector (Lx1 int)
%       model.seqlength : dna sequence length
%       model.bias : forward bias
%       model.err : error rate
%       model.nreads : number of reads generated
%
% last modified: 2012-01-20

% varargin defaults
use_input_dna = 0; % generate new dna sequence
debug = 1; % verbose output

for i = 1:length(varargin)
    if isstr(varargin{i})
        switch lower(varargin{i})
        case {'dna'} % use a previously generated dna sequence
            dna = varargin{i+1};
            use_input_dna = 1;
        case {'debug'}
        	debug = varargin{i+1};
        end
    end
end 

% step 1: generate dna
if use_input_dna
	if debug
		fprintf('using input dna sequence\n');
	end
    seqlength = length(dna);
else
	if debug
		fprintf('generating dna\n');
	end
	% generate new dna sequence
	dna = gen_dna(seqlength);
end

% generate reads
if debug
	fprintf('generating reads\n');
end
for n = 1:nreads
	% generate read sequence x and positions z 
    if debug
        fprintf('generating read %d\n',n);
    end
	[model.reads(n).x model.reads(n).z model.reads(n).e model.reads(n).t] = gen_read(dna, bias, err);
end

% model struct contains all the parameters
% simplifies passing into inference routine
model.dna = dna;
model.seqlength = seqlength;
model.bias = bias;
model.err = err;
model.nreads = nreads;
