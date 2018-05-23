function varargout = detector_cyclicCorrelation_unknownP(x,tau_0,P_max,numD)
% DETECTOR_CYCLICCORRELATION_UNKNOWNP: Obtains the detector of [1](based on
% Gardner's detector [2]) for a given set of candidate cycle frequencies.
%
%           gamma = detector_cyclicCorrelation_unknownP(x,tau_0,P_max,numD)
%           x = multivariate time series
%           tau_0 = time lag used for the covariance matrices
%           P_max = maximum expected integer part of cycle period
%           numD = grid size of fractional parts of cycle period
%
%           gamma = detector statistic
%           frac_cyclicCorr = fractional part corresponding to maximum test statistic
%
% [1] P. Urriza, E. Rebeiz, and D. Cabric, "Multiple antenna
%     cyclostationary spectrum sensing based on the cyuclic correlation
%     significance test," IEEE J. Sel. Areas Commun., vol 31, no. 11, pp.
%     2185-2195, 2013.
% [2] S.V. Schell and W.A. Gardner,"Detection of the number of
%     cyclostationary signals in unknown interference and noise," Asilomar
%     Conference on Signaals, Systems and Computers, pp.3-7, 1990.
%
gamma= zeros(P_max,1);
if nargout > 1
   frac_cyclicCorr = zeros(P_max,1); 
end
gamma_temp = zeros(numD,1);
[L,N_samples,M] = size(x);

epsilon = linspace(-0.5,0.5 - 1/numD,numD); % grid of candidate fractional parts of cycle period

for P_int=2:P_max
    P_candidate = P_int+epsilon; % candidate cycle period
    
    for ii = 1:numD
        
        % Frequency-shifted signals
        exponential = repmat(exp(-1i*2*pi*(0:N_samples-1)/P_candidate(ii)),L,1,M);
        x_shifted = x.*exponential;
        
        % Gardner/Cabric detector
        gamma_temp(ii) = detectors_Gardner(x,x_shifted,tau_0); % obtain detector statistic for each candidate cycle period
    end
    
    [gamma(P_int), id_frac] = max(gamma_temp); % choose the cycle period that maximizes the statistic 
    frac_cyclicCorr(P_int) = epsilon(id_frac);
end
varargout{1} = gamma;
if nargout > 1
   varargout{2} = frac_cyclicCorr; 
end
end