function [gamma, ln_Z] = forwback_sparse(px_z, A, At, a0, bT)
% gamma = forwback_sparse(px_z, A, At, a0, bT)
%
% Performs forward-backward message passing for HMMs using a sparse
% transition matrix (implemented in C++)
% 
% Inputs
% ------
%
%   px_z : (T x K) 
%       Observation likelihood p(x(t) | z(t)=k, theta) = px_z(t, k) 
%       given the latent state at each time point.
%
%   A : (K x K) 
%       Transition probabilities 
%
%         p(z(t+1)=l | z(t)=k, theta)  =  A(k, l)
%
%   At : (K x K) 
%       Transpose of transition probabilities
%
%   a0 : (K x 1) 
%       Boundary condition on forward sweep
%
%   bT : (K x 1) 
%       Boundary condition on backward sweep
%
% Outputs
% -------
%
%   gamma : (T x K)
%       Posterior probabilities for states
%         p(z(t)=k | x(1:T))  =  gamma(t, k)
%
%   ln_Z : float
%       Log normalization constant Z = p(x(1:T) | theta)
%
% Jan-Willem van de Meent 
% $Revision: 1.0$  $Date: 2013/02/08$

% SEE forwback_sparse.cpp for implementation
