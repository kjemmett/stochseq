function [signal pos errloc] = gen_read(dna, p, e)
% [signal pos] = gen_read(dna, p, e)
%
% Generates a random sequence, along with a set of read outs
% based on a stochastic walk from start to end.
%
%
% Inputs
% ------
%
% dna : (L x 1) int
%   Sequence to generate reads from. Each element is int from 1 to 4.
%
% p : float
%   Stepping bias, i.e. the probability of moving forward along 
%   sequence at each timestep. p=1 means deterministic forward motion,
%   p=0.5 means fully diffusive stepping.
%
% e : float
%   Error rate in the observations. When a read error occurs, the signal
%   is replaced with a random symbol
%   
%
% Outputs
% -------
%
% signal : (T x 1) int
%   Observed read-outs resulting from stochastic walk along sequence. 
%   Each element is int from 1 to 4.
%
% pos : (T x 1) int
%   Position of read at time in the sequence (int from 1 to L).
%
% errloc : (T x 1) binary
%   0 if no read error occurred
%   1 if a read error occurred
%
%
%
% $Revision: 1.02 $  $Date: 2013/02/12$

% get sequence length
L = length(dna);

%% Generate model noisy signal
signal = [];
pos = [1];
errloc = [];
t = 1;

while true
    % generate read
    if rand > e
        signal = [signal; dna(pos(t))];
        errloc = [errloc; 0];
    else
        signal = [signal; ceil(4*rand)];
        errloc = [errloc; 1];
    end
    
    % make a move
    if rand < p
        if pos(t) < L
            pos = [pos; pos(t) + 1];
        else
            % stop when seq runs beyond pos(t) > L
            break
        end
    else
        if pos(t) > 1
            pos = [pos; pos(t) - 1];
        else
            % if pos(t) drops below 1, start over
            % (this is to prevent assymmetry in the number of reads 
            % near the start of sequence compared to end)
            pos = [1];
            signal = [];
            errloc = [];
            t = 0;
        end
    end
    t = t+1;
end
