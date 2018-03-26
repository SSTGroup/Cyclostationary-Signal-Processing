function test_statistic = urriza_detector(X,tau_0,alpha_0,fs)
%URRIZA_DETECTOR detect presence of a signal in noise
% Based on P. Urriza, E. Rebeiz, and D. Cabric, "Multiple antenna 
%   cyclostationary spectrum sensing based on the cyclic correlation 
%   significance test," IEEE J. Sel. Areas Commun., vol. 31, no. 11, 
%   pp. 2185-2195, 2013.
%
%   Input:
%       X               -   input data (number of receiver antennas x number of observations)
%       tau_0           -   lag/shift (negative)
%       alpha_0         -   cyclic frequency to be tested
%       fs              -   sample frequency
%
%   Output: 
%       test_statistic  -   test statistic

N_samples = size(X,2);
M = size(X,1);

N_prime = N_samples-1-abs(tau_0); 

% obtain the two covariance matrices
R_xx_alpha_hat = bsxfun(@times,X(:,1:end+tau_0) , exp(-1j*2*pi*alpha_0*(0:N_prime)/fs) )   *ctranspose(X(:,-tau_0+1:end)) / N_prime;
R_xx_hat = X*ctranspose(X) / N_samples;

R_hat = inv(R_xx_hat)*R_xx_alpha_hat*inv(R_xx_hat)*ctranspose(R_xx_alpha_hat);

% eigenvalue decomposition
mu = eig(R_hat);

% integrated test, the eigenvalues should be real - print a warning
% otherwise
if max(abs(imag(mu))) >1e-10
    warning('Eigenvalue should not be complex')
end
mu = real(mu); % to prevent a complex test_statistic

% finally, compute the test statistic as a function of the eigenvalues
lambda = prod(1-mu);
test_statistic = -2*(N_samples-M-1)*log(lambda);

end