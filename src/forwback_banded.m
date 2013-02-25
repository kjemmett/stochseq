function  [gamma, xi, ln_Z] = forwback_banded(px_z, A, d, alpha_1, beta_T)
% [gamma, xi, ln_Z] = forwback_banded(px_z, A, d, alpha_1, beta_T)
%          
% Performs forward-backward message passing for HMMs.
% 
% Inputs
% ------
%
%   px_z : T x K 
%       Observation likelihood p(x(t) | z(t)=k, theta) = px_z(t, k) 
%       given the latent state at each time point.
%
%   A : K x L 
%       Transition probabilities for states
%
%         a(k, l) = p(z(t+1)=k+d(l) | z(t)=k, theta)
%
%       Note that a(k, l) may assign a non-zero probability to 
%       transitions outside the range 1:K. These transitions are
%       assumed to imply an invalid trajectory, and will be
%       assigned zero probability in the posteriors gamma and xi. 
%       However, the probability mass assigned to the transitions 
%       does imply a lower log likelihood p(x(1:T) | theta).
%
%   d : L x 1
%       Indices for diagonals.
%
%   alpha_1 : K x 1, optional
%       Weights for state at t=1
%
%   beta_T : K x 1, optional
%       Weights for state at t=T
%
% Outputs
% -------
%
%   gamma : T x K
%       Posterior probabilities for states
%         p(z(t)=k | x(1:T))  =  gamma(t, k)
%
%   xi : K x L
%       Posterior joint counts for states
%         xi(k, l) = sum_t p(z(t+1)=k+d(l), z(t)=k | x(1:T)) 
%
%   ln_Z : float
%       Log normalization constant Z = p(x(1:T) | theta)
%
% Jan-Willem van de Meent 
% $Revision: 1.2$  $Date: 2011/08/08$
