function [alpha_cycle,test_statistic] = cycleTest(x,P_FA,tau_list,alpha_test)
% CYCLETEST Test for cycle frequencies
% This test employs the Dandawate (1994) test statistic, but only for tau =
% 0. The test is performed using a FFT and hence the resulution of the
% cycle frequencies alpha is bound.
%
% This function has two modes: If alpha_test is set, then it evaluates the
% test statistic at that cycle frequency. If alpha_test is not set, then we
% evaluate many cycle frequency and determine which are significant.
%
%   Input:
%       x            -   signal
%       P_FA         -   desired probability of false alarm
%       tau_list     -   list of lags
%       alpha_test   -   optional: cycle frequency
%
%   Output: 
%       alpha_cycle  -   declared cycle frequencies
%       test_statistic - test statistic [18] as of Dandawate (1994)

L = 2049; % length of the window (odd)
beta = 10; % parameter of the Kaiser window
T = numel(x);

%% Step 1: C_hat
% special case tau = 0, to add more taus, we need to introduce
% minor modifications (incomplete list):
% FFT over x times shifted x (manual zero padding) -> use matrix for FFT
% covariance matrix gets larger with more taus
x_with_shifts = zeros(numel(tau_list),T);
for i = 1:numel(tau_list)
    if tau_list(i) <= 0
        x_with_shifts(i,:) = [zeros(1,abs(tau_list(i))) , x(-tau_list(i)+1:end).*conj(x(1:end+tau_list(i)))];
    elseif tau_list(i) > 0
        x_with_shifts(i,:) = [zeros(1,abs(tau_list(i))) , x(1:end-tau_list(i)).*conj(x(tau_list(i)+1:end))];
    end
end

C_hat_alpha = 1/T*fft(x_with_shifts,[],2); % columns: alphas, each row: one tau
alpha_list = 2*pi/T*linspace(0,T-1,T); % no fftshift

%% Step 2: test statistic & decision
s_list = -(L-1)/2:(L-1)/2;
W = besseli(0,beta*sqrt(1-(2*s_list/(L-1)+1-1).^2))/besseli(0,beta); % Kaiser window
W = 1/2*(W + fliplr(W)); % make W symmetric

DFT_IDX = @(n,N) mod(n-1,N)+1; % to make use of the periodicity of the DFT

if nargin == 4
    alpha_list2 = alpha_test;
else
    alpha_list2 = 1:T;
end
test_statistic = zeros(1,numel(alpha_list2));

for k_alpha = alpha_list2 % can be vectorized, at least for a single tau
    c_hat = [real(C_hat_alpha(:,k_alpha))' imag(C_hat_alpha(:,k_alpha))']; % estimate for a single alpha!
    
    weighted_C_hat_alpha = bsxfun(@times,W,C_hat_alpha(:,DFT_IDX(k_alpha + s_list,T)));
    Q = T/L*weighted_C_hat_alpha*C_hat_alpha(:,DFT_IDX(k_alpha - s_list,T)).';    
    Q_tilde = T/L*weighted_C_hat_alpha*(C_hat_alpha(:,DFT_IDX(k_alpha + s_list,T)))';    

    Sigma_hat = 1/2*[real(Q+Q_tilde) , imag(Q-Q_tilde);...
                 imag(Q+Q_tilde) , real(-Q+Q_tilde)];

    % Integrated test cases
    if norm(Sigma_hat-Sigma_hat','fro') > 1e-5
        warning('Sigma not symmetric')
    end
    if min(eig(Sigma_hat)) < -1e-5
        warning('Sigma not p.s.d.')
    end

    [P,Lambda] = eig(Sigma_hat);
    Lambda(Lambda < 1e-10) = 0; % Avoids problems for the case of zero noise
    Sigma_hat = P*Lambda*P';
    
    % Avoid inverting a singular matrix by falling back to the
    % pseudoinverse
    if rank(Sigma_hat) < size(Sigma_hat,1)
        test_statistic(k_alpha) = T*c_hat*pinv(Sigma_hat)*c_hat';
    else
        test_statistic(k_alpha) = T*c_hat/Sigma_hat*c_hat';
    end
end
test_statistic = test_statistic(alpha_test);

threshold = chi2inv(1-P_FA,2*numel(tau_list)); % N = 1 because of tau = 0 only
alpha_cycle = alpha_list(test_statistic > threshold);

end