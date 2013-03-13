function [model] = stochseq_build(seqlength, bias, err, nreads, varargin)
% [model] = stochseq_build(seqlength, bias, err, nreads, varargin)
%
% data generation for stochseq
%  1. generates random dna sequence (or uses input sequence)
%  2. builds a set of 'nreads' dna sequences
%
% use output struct as input to stochseq_infer()
%
% Arguments
% ---------
% seqlength : int
%   sequence length
% bias : scalar
%   forward bias
% err : scalar
%   error rate
% nreads : int
%   number of reads to generate
%
% Variable Arguments
% ------------------
% dna : [L 1] numeric
%   Use as input dna sequence
%   (skip gen_dna step)
% verbose : boolean
%   Print status updates (default: true)
%
% Outputs
% -------
%   model : struct containing read sequences and data parameters
%       .reads : struct containing read sequence information
%           .x : observed state sequence
%           .z : latent state sequence
%           .e : read error locations
%       .dna : [L 1] int
%           dna sequence vector
%       .seqlength : int
%           dna sequence length
%       .bias : scalar
%           forward bias
%       .err : scalar
%           error rate
%       .nreads : int
%           number of reads generated
%
% last modified: 2013-02-12

% parse varargs
ip = inputParser();
ip.StructExpand = true;
ip.addParamValue('dna', 0);
ip.addParamValue('verbose', true, @isscalar);
ip.parse(varargin{:});
args = ip.Results;

% step 1: generate dna
if args.dna==0
	if args.verbose
		fprintf('generating dna\n');
	end
	% generate new dna sequence
	dna = gen_dna(seqlength);
else
	if args.verbose
		fprintf('using input dna sequence\n');
	end
    dna = args.dna;
    seqlength = length(dna);
end

% step 2: generate reads
if args.verbose
	fprintf('generating reads\n');
end
for n = 1:nreads
	% generate read sequence x and positions z 
    if args.verbose
        fprintf('generating read %d\n',n);
    end
	[model.reads(n).x model.reads(n).z model.reads(n).e] = gen_read_triplet(dna, bias, err);
end

% model struct contains all the parameters
% simplifies passing into inference routine
model.dna = dna;
model.seqlength = seqlength;
model.bias = bias;
model.err = err;
model.nreads = nreads;
