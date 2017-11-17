function gamma = detector_Poor(x,x_shifted,alpha,Lags,window_length)
% DETECTOR_POOR: Obtains Giannakis' detector for a fixed cycle frequency
%
%
%           gamma = detector_Giannakis(x,x_shifted,alpha,L,window_length)
%           x = multivariate time series
%           x_shifted = frequency-shifted time series
%           alpha = cycle-frequency (digital domain)
%           Lags = number of lags used in the detector
%           window_length = length of the Kaiser window
%
%           gamma = detector statistic

% David Ramirez
% Signal and System Theory group
% University of Paderborn (Germany)
% 2014

[L,N] = size(x);

%% Cyclic covariance sequence
c2 = zeros(2*Lags,L);

gamma = 0;

for kk = 1:L
    c = xcorr(x_shifted(kk,:),x(kk,:),Lags-1,'biased');c = c(Lags:-1:1);
    c2(:,kk) = [real(c)' ; imag(c)'];
    
    %% Covariance matrices
    
    F_tau_omega = zeros(Lags,N);
    
    for tau = 0:Lags-1
        F_tau_omega(1+tau,:) = fft(x(kk,1:end-tau).*conj(x(kk,1+tau:end)),N);
    end
    
    Q = zeros(Lags,Lags);
    Q_ast = zeros(Lags,Lags);
    
    alpha_index = round(alpha*N);
    
    window = kaiser(window_length,10)'; window = window./sum(abs(window)); % Normalize window
    s = -(window_length-1)/2:(window_length-1)/2;
    
    for tau_m  = 0:Lags-1
        for tau_n  = 0:Lags-1
            indexes_1 = mod(alpha_index-s,N);indexes_1(indexes_1==0) = N;
            indexes_2 = mod(alpha_index+s,N);indexes_2(indexes_2==0) = N;
            Q(tau_m+1,tau_n+1) = sum(window.*F_tau_omega(1+tau_n,indexes_1).*F_tau_omega(1+tau_m,indexes_2))/N; % Due to window normalization we don't have to divide by the window_length
            Q_ast(tau_m+1,tau_n+1) = sum(window.*conj(F_tau_omega(1+tau_n,indexes_2)).*F_tau_omega(1+tau_m,indexes_2))/N;
        end
    end
    
    Sigma = [real(Q + Q_ast)/2 imag(Q - Q_ast)/2 ; imag(Q + Q_ast)/2 real(Q_ast - Q)/2];
    
    %% GLRT
    
    gamma = gamma + N*c2(:,kk)'*(Sigma\c2(:,kk));
end

