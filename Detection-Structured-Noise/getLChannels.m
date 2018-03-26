function h_channel = getLChannels(L, h_profile)
%GETLCHANNELS Sample L random channels
%   The channel model is Rayleigh fading with exponential power delay
%   profile, as given in h_profile
%
%   Input:
%       L           -   number of channels
%       h_profile   -   power delay profile (vector)
%
%   Output:
%       h_channel   -   channel taps

T = size(h_profile, 2); % number of taps
h_channel = zeros(L, T);
for l = 1:L
    h_channel(l,:) = h_profile.*randn(1,T)+1j*h_profile.*randn(1,T);
end