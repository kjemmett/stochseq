function S_new = shift_est(S,start,stop,dir)
% function S_new = shift_est(S,start,stop,dir)
%
% inputs:
%   S : sequence estimate
%   start : shift start (int)
%   stop : shift stop (int)
%   dir : direction (+1,-1)

S_new = S;

d = 2 * dir;

S_new(start+d:stop+d,:) = S(start:stop,:);

% smooth over shifted area

S_new(start-1:start+1,:) = 0.25;
S_new(stop-1:stop+1,:) = 0.25;
