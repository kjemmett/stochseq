function run_lambda_stochseq(run_id, run_path, L, p, e, N, varargin)
% function run_lambda_stochseq(run_id, run_path, L, p, e, N, varargin)

disp('=====STOCHSEQ=====');

% parse input parameters
L = str2num(L);
p = str2num(p);
e = str2num(e);
N = str2num(N);

% parse varargin
ip = inputParser();
ip.StructExpand = true;
ip.addParamValue('par', 'false');
ip.addParamValue('workers', 0);
ip.parse(varargin{:});
args = ip.Results;
args.par = str2num(args.par);

% read lambda phage sequence
disp('loading data');
[h s] = fastaread([run_path '/data/lambdaphage.fasta']);

% take only the first L bases
% cast as int array
dna = nt2int(s(1:L));

% generate model object
disp('generating model object');
model = stochseq_build(L, p, e, N, 'dna', dna, 'verbose', false);

% save model object
disp('saving model object');
savefile = [run_path '/output/' run_id '.mat'];
save(savefile, 'model');

% set up parallel workers
% if applicable
if args.par
    workers = str2num(args.workers);
    matlabpool('close', 'force', 'local')
    matlabpool('open', 'local', workers)
end

% perform inference
disp('performing inference');
inference = stochseq_infer_triplet(model, 'verbose', true);

% close parallel workers
if args.par
    matlabpool('close')
end

% save inference object
disp('saving inference object');
savefile = [run_path '/output/' run_id '.mat'];
save(savefile, 'inference');
