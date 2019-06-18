function gamma = passive_detector_cross_corr(u_s,u_r)
% PASSIVE_DETECTOR_CROSS_CORR: Cross-correlates the observations of two 
% MIMO channels as proposed in e.g. [2] and extended to multiple antennas 
% as in [1].
%
%           gamma = passive_detector_cross_corr(u_s,u_r)
%           u_s = multivariate time series observed at surveillance channel
%           u_r = multivariate time series observed at reference channel
%
%           Output: 
%           gamma = test statistic
%
[~,NMP] = size(u_s);
% Obtain sample cross-covariance matrix
S_sr = u_s*u_r'/NMP;

% Obtain the test statistic
gamma = abs(trace(S_sr'*S_sr));
