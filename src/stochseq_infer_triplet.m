function [inf_output] = stochseq_infer_triplet(model, varargin)
% function [inf_output] = stochseq_infer_triplet(model, varargin)
%
% sequence inference for stochseq
%   uses output of stochseq_build to perform
%   sequence inference using HMM model
%
% inputs:
%   model : struct containing reads, seqlength, bias, err, dna
%       see stochseq_build for full contents
%
% varagin:
%   'S0' : initial sequence estimate [default = uniform]
%   'epsilon' : likelihood convergence criterion [default = 1e-5]
%   'max_sweep' : maximium iterations of EM loop [default = 100]
%   'verbose' : verbose output? [default = true]
%
% outputs:
%   inf_output : struct containing inference ouput
%       .S : sequence estimate
%       .ll : log likelihood array
%       .??
%
% last modified: 2013-03-07

% add project path
addpath(genpath('.'));

% parse varargin
ip = inputParser();
ip.StructExpand = true;
ip.addParamValue('S0', (0.25 * ones(model.seqlength, 4)), @isnumeric); 
ip.addParamValue('epsilon', 1e-5, @isscalar);
ip.addParamValue('max_sweep', 100, @isscalar);
ip.addParamValue('verbose', true, @isscalar);
ip.parse(varargin{:});
args = ip.Results;

% pull read data from model struct
reads = model.reads;

% pull input parameters from model struct
nreads = model.nreads;
L = model.seqlength;
bias = model.bias;
dna = model.dna;
err = model.err;

% construct initial transition matrix
A = gen_transmatrix(L,bias);

% begin inference
S = args.S0;
iter = 1;
log_pX(1) = 0;
while iter <= args.max_sweep

    tic;

    if iter > 1
        log_pX(iter) = log_pX(iter - 1);
    end

	% run em for each read
	em = cell(nreads, 1);
    for n = 1:nreads
        % run em
        [em{n}.S em{n}.gamma em{n}.LpX] = em_step_triplet(reads(n).x, S, A, err);
    end
    % convert output to struct
    em = [em{:}];

	% perform hierarchical update using averaged S
	%inf_output.h(iter).em = em;
    S = h_step(em);

	% compute current iteration log likelihood
	log_pX(iter) = sum([em(:).LpX]);

	% update previous iteration log likelihood
	if (iter > 1)
        dL = abs(log_pX(iter) - log_pX(iter - 1)) / abs(log_pX(iter));
	else
		dL = NaN;
	end
	
    t = toc;
	if args.verbose
        [ignore seq] = max(S, [], 2);
        fprintf('sweep: %d  ll: %.6e  time: %.6f\n', iter, log_pX(iter), t);
        %fprintf('edit distance = %d\n',strdist(int2nt(dna'),int2nt(guess')));
	end

	% convergence check
	if (iter > 1) & (dL < args.epsilon)
        break
	end

	% proceed to next iteration
	iter = iter + 1;	
end

% diagnostics
[ignore seq_guess] = max(S, [], 2);
edit_distance = strdist(int2nt(dna'), int2nt(seq_guess'));

if args.verbose
    fprintf('converged in %d steps\n', iter);
    fprintf('log likelihood = %.6e\n', log_pX(end));
    fprintf('edit distance = %d\n', edit_distance);
end

inf_output.S = S;
inf_output.ll = log_pX;
inf_output.ed = edit_distance;
