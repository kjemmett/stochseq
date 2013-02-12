function [lambda, phi, elbo] = svb_step_sparse(x, lambda0, alpha, A, rho, N)
    % [lambda, phi, elbo] = svb_step_sparse(x, lambda0, alpha, A, rho, N)
    %
    %
    % Arguments
    % ---------
    %
    % x : [N 1] cell 
    %   Reads. Each x{n}(t) lies in range 1:4.
    %
    % lambda0 : [L D]
    %   Initial guess for Dirichlet posterior on sequence.
    %
    % alpha : [L D] numeric
    %   Dirichlet prior on sequence.
    %
    % A : [L L]
    %   Transition matrix.
    %  
    % rho : scalar
    %   Step size for stochastic update. 
    %
    % N : scalar
    %   Total number of observations in dataset
    %
    % Returns
    % -------
    %
    % lambda : [L D]
    %   Posterior on parameters. Same fields as alpha.
    %
    % phi : [N 1] cell
    %   Posterior on sequence positions. phi{n} is size [T{n} L].
    %
    % elbo : scalar
    %   Evidence lower bound
    %
    % For algorithm details see Matt Hoffman's excellent write-up: 
    % [hoffman2012] http://arxiv.org/abs/1206.7051
    %
    % (c) 2013 Jan-Willem van de Meent

    % Generative Model
    %
    % beta ~ Dir(alpha)
    % z{n}(t+1) ~ Disc(A(z{n}(t), :)) 
    % x{n}(t) | z{n}(t) ~ Disc(beta(z{n}(t),:))

    % Transition matrix A(pi)
    %
    % A = [[ 0     pi                    ]
    %      [ 1-pi  0     pi              ]
    %      [       1-pi  0     pi        ]
    %      [             ...   ...   ... ]
    %      [                   1-pi  0   ]]

    % Prior in Conjugate Exponential Form
    %
    % p(beta | alpha) 
    %   = h(beta) exp[ alpha' u(beta) - a(alpha) ]
    %
    % h(beta) = 1 ./ beta
    % u(beta) = log(beta)
    % a(alpha) = sum(log(alpha),2) - log(sum(alpha,2))

    % Probability for State Sequence in Conjugate Exponential Form
    %
    % p(x{n}, z{n} | beta) 
    %   = h(x{n}, z{n}) exp[ g(beta)' u(x{n}, z{n}) - a(g(beta)) ]
    %
    % h(x{n}, z{n}) = 1
    % g(beta) = log beta(l,d)
    % u(x{n}, z{n}) = sum_t z{n}(t,k) x{n}(t,d)
    % a(g(beta)) = 0

    % Posterior on States in Conjugate Exponential Form
    %
    % p(z{n} | x{n}, beta)
    %   =  h(z{n}) exp[eta(x{n}, beta)' u(z{n}) - a(eta(x{n}, beta))] 
    %
    % h(z{n}) = ?
    % eta(x{n}, beta) = ?
    % u(z{n}) = ?

    % get dimensions
    [L D] = size(lambda0);

    % assume batch of size one if x is not a cell array
    if not(iscell(x))
        x = {x};
    end

    % We'll need these for every read
    At = transpose(A);
    a0 = zeros(L,1);
    bT = zeros(L,1);
    a0(1) = 1;
    bT(end) = 1;

    % E_lambda [ g(beta) ] = E_lambda [ log(beta | lambda0) ]
    E_log_beta = bsxfun(@minus, psi(lambda0), psi(sum(lambda0 ,2))); 
    E_beta = exp(E_log_beta);

    % initialize sufficient statistics
    u = zeros(size(lambda0));
    for n = 1:length(x) 
        % Update variational distribution
        %
        % q(z{n} | phi) 
        %   = exp[ E_lambda [ log p(x{n}, z{n}, beta) ] ] / Z_lambda
        %   = p(x, z, exp(E_lambda [ log beta ])) / Z_lambda
        %
        % (TODO: add derivation somewhere)

        % phi(t,k) = E_lambda [ p(z{n}(t)=k | x{n}, beta) ] 
        [phi, ln_Z{n}] = forwback_sparse(E_beta(:, x{n})', A, transpose(A), a0, bT);

        % update posterior counts
        T = length(x{n});
        xind = bsxfun(@eq, x{n}(:)', (1:D)');
        u = u + sum(bsxfun(@times, ...
                           reshape(phi', [L 1 T]), ...
                           reshape(xind, [1 D T])),3);
    end

    % get number of datapoints in batch
    Nb = sum(cellfun(@length, x));

    % update lambda
    lambda = (1-rho) * lambda0 + rho * (N / Nb * u + alpha);

    % calculate lower bound
    elbo = sum([ln_Z{:}]) - sum(kl_dir(lambda, alpha),1);