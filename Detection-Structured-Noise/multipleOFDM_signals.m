function s = multipleOFDM_signals(N, P, dataOFDM, h_channel)
%MULTIPLEOFDM_SIGNALS 
%   Convenience function to sample OFDM signals for multiple antannes and
%   pass them through a channel.
%
%   Input:
%       N, P        -   number of symbols, cycle period
%       dataOFDM    -   information to be 'transmitted'
%       h_profile   -   power delay profile (vector)
%
%   Output:
%       s           -   resulting signal

L = size(h_channel, 1);
s = zeros(L, N*P);

for l = 1:L
    s(l, :) = filter(h_channel(l,:),1,dataOFDM(:,l).');
end
    