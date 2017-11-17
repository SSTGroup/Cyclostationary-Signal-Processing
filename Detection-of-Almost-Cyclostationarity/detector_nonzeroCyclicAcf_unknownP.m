function gamma = detector_nonzeroCyclicAcf_unknownP(x,Lags,window_length,P_max,numD)
% DETECTOR_NONZEROCYCLICACF_UNKNOWNP: Obtains the detector of [1](based on
% Giannakis' detector [2]) for a given set of candidate cycle frequencies.
%
%           gamma = detector_nonzeroCyclicAcf_unknownP(x,Lags,window_length,P_max,numD)
%           x = multivariate time series
%           Lags = number of lags used in the detector
%           window_length = length of the Kaiser window
%           P_max = maximum expected integer part of cycle period
%           numD = grid size of fractional parts of cycle period
%
%           gamma = detector statistic
%
% [1] J. Lunden, V. Koivunen, A. Huttunen, and H.V. Poor, "Collaborative
%     cyclostationary spectrum sensing for cognitive radio systems",  IEEE
%     Trans. Signal Process., vol. 57, no. 11, pp. 4182-4195, 2009.
% [2] A.V. Dandawate and G.B. Giannakis,"Statistical tests for presence of
%     cyclostationarity," IEEE Trans. Signal Process., vol. 42, no. 9, 1994.
%
gamma= zeros(P_max,1);
gamma_temp = zeros(numD,1);

[L,N_samples,M] = size(x);

epsilon = linspace(-0.5,0.5 - 1/numD,numD);

for P = 2:P_max
    
    P_cand = P + epsilon;
    
    % Sweep through different sampling rates D=L/M
    for ii = 1:numD
        % Frequency-shifted signals
        exponential = repmat(exp(-1i*2*pi*(0:N_samples-1)/P_cand(ii)),L,1,M);
        x_shifted = x.*exponential;
        
        % Poor/Giannakis detector
        gamma_temp(ii) = detector_Poor(x,x_shifted,1/P_cand(ii),Lags,window_length); % obtain detector statistic for each candidate cycle period
    end
    
    gamma(P)= max(gamma_temp);  % choose the cycle period that maximizes the statistic 
end
end