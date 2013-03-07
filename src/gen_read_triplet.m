function [signal pos errloc] = gen_read_triplet(dna, p, e)
% [signal pos errloc] = gen_read_triplet(dna, p, e)
%
% Performs a stochastic walk along an input dna sequence. At each
% step we see three bases: the base we are at, back one base, and
% forward one base (subject to error). As there are 4 bases, a read
% state can take 4^3 = 64 possible values. The start and end states
% will be taken to represent state 0  and 65, respectively.
%
% Inputs
% ------
%
% dna : (L x 1) int
%   Sequence from which to generate a read. Each element is int from 1 to 4.
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
% $Revision: 1.00 $  $Date: 2013/03/04$

% get sequence length
L = length(dna);

% generate model noisy signal
signal = [];
pos = [1];
errloc = [];
t = 1;

while true
    % generate read
    local_signal = [];
    local_err = [];
    if pos(t) > 1 && pos(t) < L
        for i=-1:1
            if rand > e
                local_signal = [local_signal dna(pos(t)+i)];
                local_err = [local_err 0];
            else
                local_signal = [local_signal ceil(4 * rand)];
                local_err = [local_err 1];
            end
        end
    elseif pos(t) == 1
        local_signal = [0];
        local_err = [0];
        for i=0:1
            if rand > e
                local_signal = [local_signal dna(pos(t)+i)];
                local_err = [local_err 0];
            else
                local_signal = [local_signal ceil(4 * rand)];
                local_err = [local_err 1];
            end
        end
    elseif pos(t) == L
        for i=-1:0
            if rand > e
                local_signal = [local_signal dna(pos(t)+i)];
                local_err = [local_err 0];
            else
                local_signal = [local_signal ceil(4 * rand)];
                local_err = [local_err 1];
            end
        end
        local_signal = [local_signal 0];
        local_err = [local_err 0];
    end

    signal = [signal; local_signal];
    errloc = [errloc; local_err];
    
    % make a move
    if pos(t) == L
        % stop when seq runs to end
        break
    elseif rand < p
        if pos(t) < L
            pos = [pos; pos(t) + 1];
        end
    else
        if pos(t) > 2
            pos = [pos; pos(t) - 1];
        else
            % if pos(t) drops below 2, start over
            % (this is to prevent assymmetry in the number of reads 
            % near the start of sequence compared to end)
            pos = [1];
            signal = [];
            errloc = [];
            t = 0;
        end
    end
    t = t + 1;
end
