function statistic_GLRT = detector_ACS(x,M, N_samples, P_max, numD)
% DETECTOR_ACS: Obtains the GLRT proposed in [1] for unknown cycle period.
% To this end the signal is resampled for every candidate integer part
% according to a given set of fractional parts of the cycle period such
% that the resampled signal has an integer-valued (or sufficiently small
% fractional part) cycle period. For every candidate integer period we can
% apply the detector as proposed in [1].
%
%           statistic_GLRT = detector_ACS(RxSignal,M, N_samples, P_max, numD)
%           x = multivariate time series
%           M = number of snapshots (of length N_samples)
%           N_samples = number of samples per snapshot
%           P_max = maximum expected integer part of cycle period
%           numD = grid size of fractional parts of cycle period
%
%           statistic_GLRT = detector statistic
%
% [1] D. Ramirez, P.J. Schreier, J. Via, I. Santamaria, and L.L. Scharf,
%     "Detection of multivariate cyclostationarity," IEEE Trans. Signal
%     Process., vol. 63, no. 20, p. 5395-5408, 2015
%
L = size(x,1); % Number of antennas
N = floor(N_samples*2/2.5); % Specify number of samples to consider after resampling (common for all P_int)

gamma_GLRT = zeros(numD,1);
id_GLRT = zeros(P_max-1,1);
statistic_GLRT = zeros(P_max-1,1);

epsilon = linspace(-0.5,0.5 - 1/numD,numD); % grid of candidate fractional parts of cycle period

for p = 2:P_max
    N_min = floor(N/p); % Specify number of integer mulitples of p to consider for each candidate (number of samples per snapshot = N_min*p)
    D = p./(p + epsilon); % Grid of candidate resampling rates
    
    % Sweep through different resampling rates D
    for ii = 1:numD
        [L_resam,M_resam] = rat(D(ii)); % Find upsampling and downsampling factor
        
        % Resampling
        signal_resampled = resample(x',L_resam,M_resam);
        signal_resampled = reshape(signal_resampled(1:N_min*p*M,:)',L,[],M); % Reshape long sequence into M snapshots (that are in fact non-iid realizations)
        
        % Compute likelihood function for the given sampling rate
        gamma_GLRT(ii) = detector_cyclostationarity(signal_resampled(:,1:N_min*p,:),p,M);
        clear signal_resampled
    end
    [statistic_GLRT(p),id_GLRT(p)] = min(gamma_GLRT); % Find maximum likelihood estimate of fractional part by minimizing the GLRT statistic
end