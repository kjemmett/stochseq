function [lambda, log_pX] = stochseq_infer_svb(model, varargin)
% [lambda, log_pX] = stochseq_infer_svb(model, varargin)
%
% Arguments
% ---------
% 
% model : struct
%   Input data
%   .nreads : int
%       Number of reads
%   .reads : [nreads 1] struct
%       Reads
%   .bias : scalar
%       Forward bias
%   .err : scalar
%       Error rate in reads
%   .dna : [L 1]
%       True sequence
%
% Variable Arguments
% ------------------
%
% batch : int 
%   Batch size (default: 1)
% alpha : [L 1] numeric
%   Dirichlet prior on sequence (default: ones(L, 4))
% lambda0 : [L 1] numeric
%   Initial guess for dirichlet posterior (default: homogenous)
% tau : scalar > 0
%   Discount for initial iterations (default: 1)
% kappa : scalar [0.5, 1) 
%   Forgetting factor.
% max_sweep : int
%   Maximum number of full passes over data (default: 1)
% verbose : boolean
%   Print status updates (default: true)

% variables (follows http://arxiv.org/abs/1206.7051)
% 
% lambda : dirichlet posterior on sequence
% phi : posterior on states
% beta : multinomial representation of sequence
% alpha : dirichlet prior on sequence
% tau : discount factor for early iterations (>0)
% kappa : forgetting factor [0.5, 1)

% pull data from model struct
reads = model.reads;

% count total number of datapoints
N = sum(cellfun(@length, {reads.x}));

% pull input parameters from model struct
nreads = model.nreads;
L = model.seqlength;
bias = model.bias;
dna = model.dna;
err = model.err;

% parse varargs
ip = inputParser();
ip.StructExpand = true;
ip.addParamValue('batch', 1, @isscalar);
ip.addParamValue('alpha', ones(L, 4), @isnumeric);
ip.addParamValue('lambda0', (0.25 * N + 1) * ones(L, 4), @isnumeric);
ip.addParamValue('tau', 1.0, @isscalar);
ip.addParamValue('kappa', 0.5, @isscalar);
ip.addParamValue('epsilon', 1e-5, @isscalar);
ip.addParamValue('max_sweep', 1, @isscalar);
ip.addParamValue('verbose', true, @isscalar);
ip.parse(varargin{:});
args = ip.Results;

% construct initial transition matrix
A = gen_transmatrix(L, bias);

% split up data into batches
batches = {};
n = 1;
while true
    batches{end+1} = {reads(n:min(n+args.batch, nreads)).x};
    if (n+args.batch) >= nreads
        break;
    end
    n = n + args.batch;
end

% loop for number of full sweeps over data
lambda = args.lambda0;
s = 1;
it = 1;
log_pX(1) = 0;
while s <= args.max_sweep
    if s > 1
        log_pX(s) = log_pX(s-1);
    end
	% run svb updates for each batch of reads
	for b = 1:length(batches)
        % determine step size
        rho = (it + args.tau)^(-args.kappa);
		% run em
        [lambda, phi, log_px(b)] ...
            = svb_step_sparse(batches{b}, lambda, args.alpha, A, rho, N);
        % update lower bound (todo: check this)
        log_pX(s) = (1-rho) * log_pX(s) ...
                    + rho^(it>1) * N / length(batches{b}) * log_px(b);
        % print output
        if args.verbose
            fprintf('sweep: %d  batch: %d  full lower bound: %.4e  batch lower bound: %.4e\n', ...
                    s, b, log_pX(s), log_px(b));
        end
        % increment iteration count
        it = it + 1;
	end

    % calculate relative increase of lower bound
    if (s>1)
        dL = abs(log_pX(s) - log_pX(s-1)) / abs(log_pX(s));
    else
        dL = NaN;
    end

    % print status update
	if args.verbose
        % MAP estimate for sequence S(l,d)
        [ignore seq] = max(normalize(lambda, 2), [], 2);
        % % calculate string distance w.r.t. true sequence
        % dist = strdist(int2nt(dna'),int2nt(seq'));
        fprintf('sweep: %d  rel increase: %.4e\n', s, dL);        
	end

    % convergence check
    if (s > 1) & (dL < args.epsilon)
        break
    end

    % increment sweep count
    s = s + 1;
end