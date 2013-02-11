function inf_output = stochseq_infer(model, varargin)
% function [inf_output] = stochseq_infer(model, varargin)
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
%   'epsilon' : likelihood convergence criterion [default = 1e-2]
%   'maxiter' : maximium iterations of EM loop [default = 100]
%   'debug' : verbose output? [default = 0]
%   'S0' : initial sequence estimate [default = uniform]
%   'A' : input transition matrix? [default = no]
%
% outputs:
%   inf_output : struct containing inference ouput
%     inf_output.Sequence_Inference_Distribution
%     inf_output.Sequence_Inference_Guess
%     inf_output.Sequence
%     inf_output.Edit_Distance
%     inf_output.Likelihood
%     inf_output.Sequence_Inference_Entropy
%     inf_output.Sequence_Inference_Total_Entropy
%     inf_output.Header
%     inf_output.em
%
% last modified: 2013-02-11

% pull data from model struct
reads = model.reads;

% pull input parameters from model struct
nreads = model.nreads;
L = model.seqlength;
bias = model.bias;
dna = model.dna;
err = model.err;

% varargin defaults
debug = 0; % controls console output
epsilon = 5e-3; % controls log likelihood convergence
maxiter = 100; % controls max iter of EM loop
S0 = 0.25 * ones(L,4); % controls initial sequence estimate

% parse varargin
for i = 1:length(varargin)
    if isstr(varargin{i})
        switch lower(varargin{i})
        case {'epsilon'}
            epsilon = varargin{i+1};
        case {'maxiter'}
            maxiter = varargin{i+1};
        case {'debug'}
        	debug = varargin{i+1};
        case {'s0'}
        	S0 = varargin{i+1};
    	case {'A'}
    		A = varargin{i+1};
	    end
    end
end 

% construct initial transition matrix
A = gen_transmatrix(L,bias);
%A = full(A);

% initial sequence estimate
S = S0;
inf_ent = calc_entropy(S);

% begin inference
iter = 1;

while iter < maxiter

	% run em for each read
	em = cell(nreads,1);
    for n = 1:nreads
        % run em
        [em{n}.S em{n}.gamma em{n}.LpX] = em_step_sparse(reads(n).x,S,A,err);
    end
    % convert output to struct
    em = [em{:}];

	% perform hierarchical update using averaged S
	%inf_output.h(iter).em = em;
    S = h_step(em);
    inf_output.h(iter).S = S;
	inf_ent = calc_entropy(S);
    %inf_output.h(iter).inf_ent = inf_ent;

	% update previous iteration log likelihood
	if iter>1
		sLpXprev = sLpX(iter-1);
	else
		sLpXprev = -Inf;
	end

	% compute current iteration log likelihood
	sLpX(iter) = sum([em(:).LpX]);
	
	if debug
        fprintf('\niteration %d\n',iter);
        fprintf('log likelihood = %e\n',sLpX(iter));
        [ignore guess] = max(S,[],2);
        fprintf('edit distance = %d\n',strdist(int2nt(dna'),int2nt(guess')));
        %fprintf('Dkl = %e\n',dkl(model.dna_matrix,S0));
	end

	% convergence check
	if iter > 1
		if abs(sLpX(iter) - sLpX(iter-1)) < epsilon*abs(sLpX(iter)) | sLpX(iter)==0
			break
		end		
	end

	% proceed to next iteration
	iter = iter + 1;	
end

% diagnostics
[ignore seq_guess] = max(S,[],2);
edit_distance = strdist(int2nt(dna'),int2nt(seq_guess'));
inf_ent = calc_entropy(S);
tot_ent = sum(inf_ent);

fprintf('converged in %d steps\n',iter);
fprintf('log likelihood = %d\n',sLpX(iter));
fprintf('edit distance = %d\n',edit_distance);
fprintf('total entropy = %d\n',tot_ent);

% build output struct
inf_output.S = S;
inf_output.Sequence_Inference_Distribution = S;
inf_output.Sequence_Inference_Guess = seq_guess;
inf_output.Sequence = int2nt([seq_guess]');
inf_output.Sequence_Estimate_Matrix = gen_dnamatrix(seq_guess);
inf_output.Edit_Distance = edit_distance;
inf_output.Likelihood = sLpX;
inf_output.Sequence_Inference_Entropy = inf_ent;
inf_output.Sequence_Inference_Total_Entropy = tot_ent;
inf_output.Header = [];
inf_output.em = em;
