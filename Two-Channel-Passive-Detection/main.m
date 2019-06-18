% This is a MATLAB implementation supporting the papers
%
% "Two-Channel Passive Detection Exploiting Cyclostationarity"
% by Stefanie Horstmann, David Ramirez and Peter J. Schreier,
% accepted for EUSIPCO 2019
%
% "Two-Channel Passive Detection of Cyclostationary Signals"
% by Stefanie Horstmann, David Ramirez and Peter J. Schreier,
% submitted to IEEE Transaction on Signal Processing, June 2019
%
% ## ----------------------------------------------------------------------------
% ##
% ##   File: main.m
% ##   Copyright (c) <2019> Signal and System Theory Group,
% ##                        Univ. of Paderborn, http://sst.upb.de
% ##                        https://github.com/SSTGroup/Cyclostationary-Signal-Processing
% ##
% ##   Permission is hereby granted, free of charge, to any person
% ##   obtaining a copy of this software and associated documentation
% ##   files (the "Software"), to deal in the Software without restriction,
% ##   including without limitation the rights to use, copy, modify and
% ##   merge the Software, subject to the following conditions:
% ##
% ##   1.) The Software is used for non-commercial research and
% ##       education purposes.
% ##
% ##   2.) The above copyright notice and this permission notice shall be
% ##       included in all copies or substantial portions of the Software.
% ##
% ##   3.) Publication, distribution, sublicensing, and/or selling of
% ##       copies or parts of the Software requires special agreements
% ##       with the Signal and System Theory Group, University of Paderborn,
% ##       and is in general not permitted.
% ##
% ##   4.) Modifications or contributions to the Software must be
% ##       published under this license.
% ##
% ##   5.) Any publication that was created with the help of this Software
% ##       or parts of this Software must include a citation of the paper
% ##       referenced above.
% ##
% ##   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
% ##   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
% ##   OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
% ##   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
% ##   HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
% ##   WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
% ##   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
% ##   OTHER DEALINGS IN THE SOFTWARE.
% ##
% ##   Persons using the Software are encouraged to notify the
% ##   Signal and System Theory Group at the University of Paderborn
% ##   about bugs.
% ##
% ##
% ##   Author: Stefanie Horstmann, <stefanie.horstmann@sst.upb>
% ##   Date: 14.6.2019
% ##   Version: v1
% ##
% ##   At last compiled with Matlab 2018b
% ##   Required Toolboxes:
% ##        Communications System Toolbox
% ##        Signal Processing Toolbox
% ##        Statistics Toolbox
% ## ----------------------------------------------------------------------------
%% References
% [1] I. Santamaria, L. L. Scharf, J. Via, H. Wang, and Y. Wang, “Passive
%     detection of correlated subspace signals in two MIMO channels,” IEEE
%     Trans. on Signal Process., vol. 65, no. 20, pp. 5266–5280, Oct. 2017.
% [2] J. Liu, H. Li, and B. Himed, “On the performance of the 
%     cross-correlation detector for passive radar applications,” 
%     Signal Process., vol. 113, pp. 32–37, 2015.
%

clear;close all;clc;
RandStream.setGlobalStream(RandStream('mcg16807','seed',sum(100*clock)))

%% Parameters
P = 4;                  % cycle period

L = 2;                  % Number of receive antennas at surveillance channel
L_I = 2;                 % Number of transmit antennas
fs = 1.2e6;              % Sampling frequency
N = 128;                  % Number of symbols
T = N*P;
M = 32;                  % Number of snapshots

Tc = 10;                 % Number of channel taps at sampling rate
alpha_channel = 1e5;     % Parameter of the exponential power delay profile
Tn = 20;                 % Order of the MA noise at the sampling rate

is_RC = 1;               % 1: Raised cosine shaping filters with roll-off factor = 1
                         % 0: Rectangular shaping filters

nsim = 10000;             % Number of Monte Carlo simulations

% Specify the desired SNR (here we choose SNR = SNR_s = SNR_r)
SNR = -20:-5;    % Signal to noise ratio

% Fix desired probability of false alarm
pfa = 0.01;

% Specify SNR for which we want to plot the ROC curve
SNR_plot = -10;

%% Other parameters
% Specify shaping filter
if is_RC == 1
    % RC filter
    StopBandAtn = 60;RollOff = 1;
    Hd = design(fdesign.pulseshaping(P,'Raised Cosine','Ast,Beta',StopBandAtn,RollOff));
    h_shaping = Hd.numerator;h_shaping = h_shaping/norm(h_shaping)*sqrt(P);
else
    % Rectangular filter
    h_shaping = ones(1,P);
end

% Frequency-selective channel
Delay_Paths = (0:Tc*P-1)/fs;
pdp = exp(-alpha_channel*Delay_Paths);pdp = pdp/sum(pdp);

H = comm.RayleighChannel('SampleRate',fs,'MaximumDopplerShift',0,'PathDelays',Delay_Paths,'AveragePathGains',10*log10(pdp));

N_samples = T*M;    % Number of samples

sigma_squared = 10.^(-SNR/10);  % Noise variance

%% -------- Obtain realizations for different SNR (SNR_s = SNR_r) ------ %%
gamma_H1_GLRT = zeros(length(SNR),nsim);
gamma_H0_GLRT = zeros(length(SNR),nsim);
gamma_H1_LMPIT_LS = zeros(length(SNR),nsim);
gamma_H0_LMPIT_LS = zeros(length(SNR),nsim);
gamma_H1_LMPIT_LSR = zeros(length(SNR),nsim);
gamma_H0_LMPIT_LSR = zeros(length(SNR),nsim);
gamma_H1_subspace = zeros(length(SNR),nsim);
gamma_H0_subspace = zeros(length(SNR),nsim);
gamma_H1_cross_corr = zeros(length(SNR),nsim);
gamma_H0_cross_corr = zeros(length(SNR),nsim);

for kkk = 1:length(SNR)
    for kk = 1:nsim
        
        if rem(kk-1,ceil(nsim/5)) == 0
            disp(['SNR = ' num2str(SNR(kkk)) ', Simulation ' num2str(kk) ' of ' num2str(nsim)])
        end
        
        % ---- Transmitted signal -----
        
        bits = rand(L_I,2*N*M)>0.5; % Bit generation
        symbols = (2*bits(:,1:2:end) - 1 + 1i*(2*bits(:,2:2:end) - 1))/sqrt(2);  % Mapping
        
        % Shaping
        s_temp = zeros(L_I,N*P*M);
        s_temp(:,1:P:end) = symbols;
        s_temp = filter(h_shaping,1,s_temp.').';  % Shaping: Delete transient
        s = reshape(s_temp,L_I,N*P*M);
        
        % Channel effect. We assume a square channel
        
        s_tx_s = zeros(L,T*M);
        s_tx_r = zeros(L,T*M);
        for nr = 1:L
            for nt = 1:L_I
                % surveillance channel
                s_tx_s(nr,:) = s_tx_s(nr,:) + H(s(nt,1:N_samples)')';
                reset(H)
                
                % reference channel
                s_tx_r(nr,:) = s_tx_r(nr,:) + H(s(nt,1:N_samples)')';
                reset(H)
            end
        end
        
        % ---- Received signal -----
        
        % Colored Gaussian noise, correlated among antennas,
        % independent between reference and surveillance array.
        
        % Coefficients for moving average filter
        noise_coeffs_s = (randn(L,Tn) + 1i*randn(L,Tn))/sqrt(2*Tn);
        noise_coeffs_r = (randn(L,Tn) + 1i*randn(L,Tn))/sqrt(2*Tn);
        
        noise_s = (randn(L,N_samples) + 1i*randn(L,N_samples)); % noise at surveillance channel
        noise_r = (randn(L,N_samples) + 1i*randn(L,N_samples)); % noise at reference channel
        
        for n = 1:L
            noise_s(n,:) = filter(noise_coeffs_s(n,:),1,noise_s(n,:));
            noise_r(n,:) = filter(noise_coeffs_r(n,:),1,noise_r(n,:));
        end
        
        sigmasqrt_s = randn(L) + 1j*randn(L);
        sigmasqrt_r = randn(L) + 1j*randn(L);
        
        % Scale observations to ensure desired SNR
        trace_Rs = trace(s_tx_s*s_tx_s'/N_samples);
        noise_s = sigmasqrt_s*noise_s/sqrt(trace((sigmasqrt_s*noise_s)*(sigmasqrt_s*noise_s)'/N_samples));
        noise_s = sqrt(trace_Rs*sigma_squared(kkk))*noise_s;
        
        trace_Rr = trace(s_tx_r*s_tx_r'/N_samples);
        noise_r = sigmasqrt_r*noise_r/sqrt(trace((sigmasqrt_r*noise_r)*(sigmasqrt_r*noise_r)'/N_samples));
        noise_r = sqrt(trace_Rr*sigma_squared(kkk))*noise_r;
        
        u_s = s_tx_s + noise_s; % signal plus noise at surveillance channel
        u_r = s_tx_r + noise_r; % signal plus noise at reference channel
        
        %% Compute the test statistics
        % GLRT and LMPIT-inspired test
        [gamma_H1_GLRT(kkk,kk),gamma_H1_LMPIT_LSR(kkk,kk),~]...
            = passive_detector_cyclostationarity(u_s,u_r,P,M);
        [gamma_H0_GLRT(kkk,kk),gamma_H0_LMPIT_LSR(kkk,kk,:),~]...
            = passive_detector_cyclostationarity(noise_s,u_r,P,M);
        
        % Correlated subspace detector of [1]
        [gamma_H1_subspace(kkk,kk,:)] ...
            = passive_detector_corr_subspace(u_s,u_r,L_I);
        [gamma_H0_subspace(kkk,kk,:)] ...
            = passive_detector_corr_subspace(noise_s,u_r,L_I);
        
        % Cross-correlation detector of e.g. [2]
        gamma_H1_cross_corr(kkk,kk,:) ...
            = passive_detector_cross_corr(u_s,u_r);
        gamma_H0_cross_corr(kkk,kk,:) ...
            = passive_detector_cross_corr(noise_s,u_r);
        
    end
end

%% -------- Evaluate performance -------- %%

% -------Probability of detection vs. SNR -----------

P_FA_GLRT = zeros(length(SNR),1);
P_D_GLRT = zeros(length(SNR),1);
P_FA_LMPIT_inspired = zeros(length(SNR),1);
P_D_LMPIT_LSR = zeros(length(SNR),1);
P_FA_subspace = zeros(length(SNR),1);
P_D_subspace = zeros(length(SNR),1);
P_FA_crosscorr = zeros(length(SNR),1);
P_D_cross_corr = zeros(length(SNR),1);

% Estimate probabiliy of detection for given probability of false alarm
for kkk = 1:length(SNR)
    [f_H0_GLRT(kkk,:),x_H0_GLRT(kkk,:)] = ecdf(reshape(gamma_H0_GLRT(kkk,:),[],1));
    [f_H0_LMPIT_LSR(kkk,:),x_H0_LMPIT_LSR(kkk,:)] = ecdf(reshape(gamma_H0_LMPIT_LSR(kkk,:),[],1));
    [f_H0_subspace(kkk,:),x_H0_subspace(kkk,:)] = ecdf(reshape(gamma_H0_subspace(kkk,:),[],1));
    [f_H0_cross_corr(kkk,:),x_H0_cross_corr(kkk,:)] = ecdf(reshape(gamma_H0_cross_corr(kkk,:),[],1));
    
    idx_pf_GLRT = find(f_H0_GLRT(kkk,:) <= pfa,1,'last'); thres_GLRT = x_H0_GLRT(kkk,idx_pf_GLRT);
    idx_pf_LMPIT_LSR = find(f_H0_LMPIT_LSR(kkk,:) >= 1-pfa,1,'first'); thres_LMPIT_LSR = x_H0_LMPIT_LSR(kkk,idx_pf_LMPIT_LSR);
    idx_pf_subspace = find(f_H0_subspace(kkk,:) >= 1-pfa,1,'first'); thres_subspace = x_H0_subspace(kkk,idx_pf_subspace);
    idx_pf_cross_corr = find(f_H0_cross_corr(kkk,:) >= 1-pfa,1,'first'); thres_cross_corr = x_H0_cross_corr(kkk,idx_pf_cross_corr);
    
    P_D_GLRT(kkk) = sum(gamma_H1_GLRT(kkk,:) <= thres_GLRT)/nsim;
    P_D_LMPIT_LSR(kkk) = sum(gamma_H1_LMPIT_LSR(kkk,:) >= thres_LMPIT_LSR)/nsim;
    P_D_subspace(kkk) = sum(gamma_H1_subspace(kkk,:) >= thres_subspace)/nsim;
    P_D_cross_corr(kkk) = sum(gamma_H1_cross_corr(kkk,:) >= thres_cross_corr)/nsim;
    
end
% Plot probability of detection vs. SNR
figure
hold on
semilogy(SNR,P_D_GLRT)
semilogy(SNR,P_D_LMPIT_LSR)
semilogy(SNR,P_D_subspace)
semilogy(SNR,P_D_cross_corr)
legend('GLRT','LMPIT-inspired','Corr. sub.', 'Cross corr.')
xlabel('SNR')
ylabel('p_d')
ylim([0,1])


% ---------- Obtain ROC curve for a particular SNR -------------
id_SNR = find(SNR == SNR_plot);

P_D_GLRT_ROC = zeros(1,nsim+1);
P_D_LMPIT_LSR_ROC = zeros(1,nsim+1);
P_D_subspace_ROC = zeros(1,nsim+1);
P_D_cross_corr_ROC = zeros(1,nsim+1);
for kk = 1:nsim+1
    P_D_GLRT_ROC(kk) = length(find(gamma_H1_GLRT(id_SNR,:) <= x_H0_GLRT(id_SNR,kk)))/nsim;
    P_D_LMPIT_LSR_ROC(kk) = length(find(gamma_H1_LMPIT_LSR(id_SNR,:) > fliplr(x_H0_LMPIT_LSR(id_SNR,kk))))/nsim;
    P_D_subspace_ROC(kk) = length(find(gamma_H1_subspace(id_SNR,:) > fliplr(x_H0_subspace(id_SNR,kk))))/nsim;
    P_D_cross_corr_ROC(kk) = length(find(gamma_H1_cross_corr(id_SNR,:) > fliplr(x_H0_cross_corr(id_SNR,kk))))/nsim;
end
% Plot ROC curve
figure
plot(f_H0_GLRT(id_SNR,:),P_D_GLRT_ROC)
hold on
plot(1-f_H0_LMPIT_LSR(id_SNR,:),P_D_LMPIT_LSR_ROC)
plot(1-f_H0_subspace(id_SNR,:),P_D_subspace_ROC)
plot(1-f_H0_cross_corr(id_SNR,:),P_D_cross_corr_ROC)
legend('GLRT','LMPIT-inspired','Corr. sub.', 'Cross corr.')
xlabel('p_{fa}')
ylabel('p_d')
ylim([0,1])
