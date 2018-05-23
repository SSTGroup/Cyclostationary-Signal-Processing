% This is a MATLAB implementation supporting the papers
%
% "Detection of Almost-Cyclostationarity: An Approach Based on a Multiple 
% Hypothesis Test" by Stefanie Horstmann, David Ramirez and 
% Peter J. Schreier, Proc. Asilomar Conf. Signals, Systems and
% Computers, November 2017
%
% "Joint Detection of Almost-Cyclostationary Signals and Estimation of 
% their Cycle Period" by Stefanie Horstmann, David Ramirez and 
% Peter J. Schreier, Submitted in IEEE Signal Processing Letters, May 2017
% 
% ## ----------------------------------------------------------------------------
% ##
% ##   File: main.m 
% ##   Copyright (c) <2017> Signal and System Theory Group, 
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
% ##   Date: 20.5.2018 
% ##   Version: v2  
% ##
% ##   At last compiled with Matlab 2017b
% ##   Required Toolboxes:
% ##        Communications System Toolbox
% ##        Signal Processing Toolbox
% ##        Statistics Toolbox
% ## ----------------------------------------------------------------------------
%% References
% [1] D. Ramirez, P.J. Schreier, J. Via, I. Santamaria, and L.L. Scharf,
%     "Detection of multivariate cyclostationarity," IEEE Trans. Signal
%     Process., vol. 63, no. 20, p. 5395-5408, 2015
% [2] P. Urriza, E. Rebeiz, and D. Cabric, "Multiple antenna
%     cyclostationary spectrum sensing based on the cyclic correlation
%     significance test," IEEE J. Sel. Areas Commun., vol 31, no. 11, pp.
%     xs2185-2195, 2013.
% [3] J. Lunden, V. Koivunen, A. Huttunen, and H.V. Poor, "Collaborative
%     cyclostationary spectrum sensing for cognitive radio systems",  IEEE
%     Trans. Signal Process., vol. 57, no. 11, pp. 4182-4195, 2009.
%

clear;close all;clc;
RandStream.setGlobalStream(RandStream('mcg16807','seed',sum(100*clock)))

%% Parameters
P = 4.2;                 % Cycle period
P_max = 10;              % Maximum expected integer cycle period

numD = 20;               % Size of grid of fractional parts

L = 2;                   % Number of antennas
N = 100;                 % Number of symbols
M = 25;                  % Number of snapshots

fs = 1.2e6;              % Sampling frequency
Tc = 10;                 % Number of channel taps at sampling rate
alpha_channel = 1e5;     % Parameter of the exponential power delay profile
Tn = 20;                 % Order of the MA noise at the sampling rate

is_RC = 1;               % 1: Raised cosine shaping filters. 0: Rectangular shaping filters

nsim = 400;             % Number of Monte Carlo simulations 

SNR = -8:-2;           % Signal to noise ratio

% Competitors
Lags_nonzeroCyclicAcf = 4;  % Number of lags for non-zero cycic ACF detector
window_length = 1025;    % Length of the Kaiser window

alpha_threshold = 5e-02; % probability of false alarm (plot probability of detection vs. SNR)

bta = 5; % Window parameter of kaiser window for interpolation filter
N_len = 10; % Length of kaiser window for the interpolation filter
%% Other parameters
[P_gen, P_sam] = rat(P); % P_gen is the cycle period of the CS process that
% is subsampled by P_sam to obtain the desired ACS process with cycle 
% period P

P_int = round(P); % integer part of the cycle period
N_samples = floor(P)*N*M; % Total number of samples
Rs = fs/P;          % Symbol rate

% Choose an appropriate lag tau_0 for the cyclic correlation detector
if is_RC == 1
    % Raised cosine filter
    tau_0 = 0;
else
    % Rectangular filter
    tau_0 = round(P/2);
end

if is_RC == 1
    % Raised-cosine filter
    StopBandAtn = 60;RollOff = 1;    
    Hd = design(fdesign.pulseshaping(P_gen,'Raised Cosine','Ast,Beta',StopBandAtn,RollOff));
    h_shaping = Hd.numerator;h_shaping = h_shaping/norm(h_shaping)*sqrt(P_gen);
else
    % Rectangular filter
    h_shaping = ones(1,P_gen);
end


% Noise variance
sigma2 = 10.^(-SNR/10);

% Exponential for frequency shift of the competing techniques

exponential = repmat(exp(-1i*2*pi*Rs*(0:N_samples-1)/fs),L,1);

% Frequency-selective channel

Delay_Paths = (0:Tc*P_int-1)/fs;
pdp = exp(-alpha_channel*Delay_Paths);pdp = pdp/sum(pdp);

average = sum(Delay_Paths.*pdp);
rms = sqrt(sum((Delay_Paths-average).^2.*pdp));

H = rayleighchan(1/fs,0,Delay_Paths,10*log10(pdp));
H.ResetBeforeFiltering = 1;

% Obtain parameters for filter bank
% Investigate resampling rates
epsilon = linspace(-0.5,0.5 - 1/numD,numD)';
p_vec = (2:P_max); p_vec = repmat(p_vec,numD,1); p_vec = reshape(p_vec,[],1);
p_frac_vec = p_vec + repmat(epsilon,P_max-1,1);
resam_rates_vec = p_vec./(p_frac_vec);
resam_rates_vec(resam_rates_vec == 1) = [];

[L_resam,M_resam] = rat(resam_rates_vec, 1e-12);
[L_unique,id_L_resam,ic_L] = unique(L_resam);

factor_matrix = zeros(length(id_L_resam),1);
for ii = 1:length(id_L_resam)
    factor_temp =  sort(factor(L_unique(ii)),'descend');
    if length(factor_temp) > size(factor_matrix,2)
        factor_matrix = [factor_matrix,zeros(length(id_L_resam),abs(length(factor_temp) - size(factor_matrix,2)))];
    else
        factor_temp = [factor_temp,zeros(1,abs(length(factor_temp) - size(factor_matrix,2)))];
    end
    factor_matrix(ii,:) = factor_temp;
end

%% -------- Obtain realizations for different SNR -------- %%
statistic_GLRT_H0 = zeros(length(SNR),nsim,P_max);
statistic_GLRT_H1 = zeros(length(SNR),nsim,P_max);
frac_GLRT = zeros(length(SNR),nsim,P_max);

statistic_cyclicCorrelation_H0 = zeros(length(SNR),nsim,P_max);
statistic_cyclicCorrelation_H1 = zeros(length(SNR),nsim,P_max);
frac_cyclicCorrelation = zeros(length(SNR),nsim,P_max);

statistic_nonzeroCyclicAcf_H0 = zeros(length(SNR),nsim,P_max);
statistic_nonzeroCyclicAcf_H1 = zeros(length(SNR),nsim,P_max);
frac_nonzeroCyclicAcf = zeros(length(SNR),nsim,P_max);

for kkk = 1:length(SNR)

    for kk = 1:nsim
                
        if rem(kk-1,ceil(nsim/5)) == 0
            disp(['SNR = ' num2str(SNR(kkk)) ', Simulation ' num2str(kk) ' of ' num2str(nsim)])
        end
        
        % Transmitted signal
        bits = rand(L,2*(N+2)*M)>0.5; % Bit generation
        symbols = (2*bits(:,1:2:end) - 1 + 1i*(2*bits(:,2:2:end) - 1))/sqrt(2);  % Mapping (QPSK)
        
        % In order to generate an almost-cyclostationary signal we shape
        % the symbol sequence with a filter of length P_gen. This sequence
        % is sub-sampled afterwards with P_sam in order to obtain the
        % desired cycle period P = P_gen/P_sam.
        
        % Shaping
        s_temp = zeros(L,(N+2)*P_gen*M);
        s_temp(:,1:P_gen:end) = symbols;
        s_temp = filter(h_shaping,1,s_temp.').';
        s_temp = reshape(s_temp,L,(N+2)*P_gen*M);    
        
        % Sample signal with non-submultiple
        s = zeros(L,N_samples);
        for l = 1:L
            random_start = randi(P_gen); % non-synchronized sequence
            temp = s_temp(l,random_start:P_sam:end-(P_gen-random_start));
            s(l,:) = temp(1,1:N_samples);
        end
        
        % Channel effect. We assume a square channel
        s_tx = zeros(L,N_samples);
        for nr = 1:L
            for nt = 1:L
                s_tx(nr,:) = s_tx(nr,:) + filter(H,s(nt,:));
            end
        end
            
        % Received signal
        noise_coeffs = (randn(L,Tn) + 1i*randn(L,Tn))/sqrt(2*Tn); % MA filter coefficients
        
        x_H0 = sqrt(sigma2(kkk)/2)*(randn(L,N_samples) + 1i*randn(L,N_samples)); 
        for nr = 1:L
            x_H0(nr,:) = filter(noise_coeffs(nr,:),1,x_H0(nr,:)); % H0: noise only
        end
        x_H1 = s_tx + x_H0; % H1: signal plus noise
        
        % Proposed detector based on the GLRT of [1]
%         statistic_GLRT_H1(kkk,kk,:) = detector_ACS(x_H1,M, N_samples/M, P_max, numD);
%         statistic_GLRT_H0(kkk,kk,:) = detector_ACS(x_H0,M, N_samples/M, P_max, numD);
        [statistic_GLRT_H1(kkk,kk,:),~, frac_GLRT(kkk,kk,:)] = detector_ACS_filterbank(x_H1,M, N_samples/M, P_max, numD, factor_matrix, N_len, bta);
        statistic_GLRT_H0(kkk,kk,:) = detector_ACS_filterbank(x_H0,M, N_samples/M, P_max, numD, factor_matrix, N_len, bta);
        
        % Competing detector based on correlation of the signal itself and
        % its frequency-shifted version [2]
        [statistic_cyclicCorrelation_H1(kkk,kk,:), frac_cyclicCorrelation(kkk,kk,:)] = detector_cyclicCorrelation_unknownP(x_H1,tau_0,P_max,numD);
        statistic_cyclicCorrelation_H0(kkk,kk,:) = detector_cyclicCorrelation_unknownP(x_H0,tau_0,P_max,numD);
        
        % Competing detector based on non-zero cyclic autocorrelation
        % function [3]
        [statistic_nonzeroCyclicAcf_H1(kkk,kk,:), frac_nonzeroCyclicAcf(kkk,kk,:)] = detector_nonzeroCyclicAcf_unknownP(x_H1,Lags_nonzeroCyclicAcf,window_length,P_max,numD);
        statistic_nonzeroCyclicAcf_H0(kkk,kk,:) = detector_nonzeroCyclicAcf_unknownP(x_H0,Lags_nonzeroCyclicAcf,window_length,P_max,numD);
                
    end
end

%% -------- Evaluate performance -------- %%
% Other parameters
N_samples = floor(P)*N;
T_hat = floor(N_samples*2/2.5);

% Thresholds of competing techniques
DoF_cyclicCorr = 2*L^2; % degrees of freedom of cyclic correlation detector
threshold_cyclicCorr = chi2inv((1-alpha_threshold)^(1/(numD*(P_max-1))),DoF_cyclicCorr); % numD*(P_max-1)-th order statistic

DoF_nonzeroCyclicAcf = 2*Lags_nonzeroCyclicAcf*L; % degrees of freedom of non-zero cyclic ACF detector
threshold_nonzeroCyclicAcf = chi2inv((1-alpha_threshold)^(1/(numD*(P_max-1))),DoF_nonzeroCyclicAcf); % numD*(P_max-1)-th order statistic

P_FA_GLRT = zeros(length(SNR),1);
P_D_GLRT = zeros(length(SNR),1);
P_estimate = zeros(length(SNR),nsim);
P_D_detection_estimation_GLRT = zeros(length(SNR),1);
P_FA_cyclicCorr = zeros(length(SNR),1);
P_D_cyclicCorr = zeros(length(SNR),1);
P_estimate_cyclicCorr = zeros(length(SNR),nsim);
P_D_detection_estimation_cyclicCorr = zeros(length(SNR),1);
P_FA_nonzeroCyclicAcf = zeros(length(SNR),1);
P_D_nonzeroCyclicAcf = zeros(length(SNR),1);
P_estimate_nonzeroCyclicAcf = zeros(length(SNR),nsim);
P_D_detection_estimation_nonzeroCyclicAcf = zeros(length(SNR),1);
    
for kkk = 1: length(SNR)
    % Estimate p-values for each candidate integer part P_int
    statistic_GLRT_1 = reshape(statistic_GLRT_H0(kkk,:,1:P_max),nsim,P_max);
    statistic_GLRT_2 = reshape(statistic_GLRT_H1(kkk,:,1:P_max),nsim,P_max);
    
    p_value_GLRT_H0 = zeros(nsim,P_max-1);
    p_value_GLRT_H1 = zeros(nsim,P_max-1);
    dist_GLRT_H0 = zeros(nsim,P_max-1);
    dist_GLRT_H1 = zeros(nsim,P_max-1);
    
    % It is recommended to draw the samples from the product of beta random variables off-line. 
    for p = 2:P_max
        % Draw samples from beta distribution and form the particular sums 
        % of the logarithms
        N_temp = floor(T_hat/p);
        i_prod = p;
        n_prod = L;
        beta = repmat(((2:i_prod)'-1)*n_prod,1,n_prod);
        alpha = M-beta - repmat(0:n_prod-1,i_prod-1,1);
        logY= zeros(10000,1); % the more samples you draw from the beta 
        % distribution the more accurate the estimates of your p-values get
        for nn = 1:10000
            for mm = 1:N_temp
                Y = betarnd(alpha,beta);
                logY(nn,mm) = -sum(sum(log(Y),2),1);
            end
        end
        logY_G1 = sum(logY,2);
        
        % Estimate 1st-order statistic from sum of log(beta)
        [cdf_beta,x_beta] = ecdf(-logY_G1);
        cdf_OS = 1-(1-cdf_beta).^numD;
                
        for k = 1:nsim
            % Estimate p-value for GLRT
            [~,id_1] = min(abs(x_beta-statistic_GLRT_1(k,p)));
            [~,id_2] = min(abs(x_beta-statistic_GLRT_2(k,p)));
            
            p_value_GLRT_H0(k,p-1) = cdf_OS(id_1);
            p_value_GLRT_H1(k,p-1) = cdf_OS(id_2);
        end
    end
    
    %% Multiple Hypothesis Test
    % Bonferroni/Holm multiple test procedure
    threshold = alpha_threshold/(P_max-1); 
    min_p_value_H0 = min(p_value_GLRT_H0,[],2);
    [min_p_value_H1, id_min] = min(p_value_GLRT_H1,[],2);
    P_FA_GLRT(kkk) = sum(min_p_value_H0 < threshold)/nsim;
    P_D_GLRT(kkk) = sum(min_p_value_H1 < threshold)/nsim;
    
    P_int_GLRT = id_min + 1;
    P_frac_GLRT = zeros(nsim,1);
    for ii = 1:nsim
        P_frac_GLRT(ii) = frac_GLRT(kkk,ii,id_min(ii)+1);
    end
    P_estimate(kkk,:) = P_int_GLRT + P_frac_GLRT;
    % Estimate probablity of jointly detecting ACS and correctly estimating
    % the cycle period
    P_D_detection_estimation_GLRT(kkk) = sum(((P_estimate(kkk,:) >= P - 1/numD) & (P_estimate(kkk,:) <= P + 1/numD))' & (min_p_value_H1 < threshold))/nsim;
    
    %% Competing techniques
    % Obtain the maximum (with respect to the grid of possible cycle period) test statistic of the cyclic correlation detector
    statistic_H0_cyclicCorr = max(2*M*N_samples*statistic_cyclicCorrelation_H0(kkk,:,1:P_max),[],3);
    [statistic_H1_cyclicCorr, P_int_cyclicCorr] = max(2*M*N_samples*statistic_cyclicCorrelation_H1(kkk,:,1:P_max),[],3);
    
    % Estimate the fractional part of the cycle period
    P_frac_cyclicCorr = zeros(nsim,1);
    for ii = 1:nsim
        P_frac_cyclicCorr(ii) = frac_cyclicCorrelation(kkk,ii,P_int_cyclicCorr(ii));
    end
    
    P_D_cyclicCorr(kkk) = length(find(statistic_H1_cyclicCorr > threshold_cyclicCorr))/nsim;
    P_FA_cyclicCorr(kkk) = length(find(statistic_H0_cyclicCorr > threshold_cyclicCorr))/nsim;
    
    P_estimate_cyclicCorr(kkk,:) = P_int_cyclicCorr' + P_frac_cyclicCorr;
    % Estimate probablity of jointly detecting ACS and correctly estimating
    % the cycle period
    P_D_detection_estimation_cyclicCorr(kkk) = sum(((P_estimate_cyclicCorr(kkk,:) >= P - 1/numD) & (P_estimate_cyclicCorr(kkk,:) <= P + 1/numD))...
        & (statistic_H1_cyclicCorr > threshold_cyclicCorr))/nsim;
    
    % Obtain the maximum (with respect to the grid of possible cycle period) test statistic of the non-zero cyclic ACF detector
    statistic_H0_nonzeroCyclicAcf = max(statistic_nonzeroCyclicAcf_H0(kkk,:,1:P_max),[],3);
    [statistic_H1_nonzeroCyclicAcf, P_int_nonzeroCyclicAcf] = max(statistic_nonzeroCyclicAcf_H1(kkk,:,1:P_max),[],3);
    
    % Estimate the fractional part of the cycle period
    P_frac_nonzeroCyclicAcf = zeros(nsim,1);
    for ii = 1:nsim
        P_frac_nonzeroCyclicAcf(ii) = frac_nonzeroCyclicAcf(kkk,ii,P_int_nonzeroCyclicAcf(ii));
    end
    
    P_D_nonzeroCyclicAcf(kkk) = length(find(statistic_H1_nonzeroCyclicAcf > threshold_nonzeroCyclicAcf))/nsim;
    P_FA_nonzeroCyclicAcf(kkk) = length(find(statistic_H0_nonzeroCyclicAcf > threshold_nonzeroCyclicAcf))/nsim;
    
    P_estimate_nonzeroCyclicAcf(kkk,:) = P_int_nonzeroCyclicAcf' + P_frac_nonzeroCyclicAcf;
    % Estimate probablity of jointly detecting ACS and correctly estimating
    % the cycle period
    P_D_detection_estimation_nonzeroCyclicAcf(kkk) = sum(((P_estimate_nonzeroCyclicAcf(kkk,:) >= P - 1/numD) & (P_estimate_nonzeroCyclicAcf(kkk,:) <= P + 1/numD))...
        & (statistic_H1_nonzeroCyclicAcf > threshold_nonzeroCyclicAcf))/nsim;
    
end
%% Plot probability of missed detection vs. SNR
figure
hold on
semilogy(SNR(1:end),1-P_D_GLRT(1:end),'b-o')
semilogy(SNR(1:end),1-P_D_nonzeroCyclicAcf(1:end),'g-*')
semilogy(SNR(1:end),1-P_D_cyclicCorr(1:end),'m-+')
legend('Proposed detector','Based on [3]','Based on [2]')
xlabel('SNR')
ylabel('p_m')
ylim([0,1])
title(['P=',num2str(P),', N=',num2str(N),', L=',num2str(L), ', M=',...
    num2str(M), ', D=',num2str(numD)])

%% Plot probability of jointly detecting and correctly estimating the cycle period vs. SNR
figure
hold on
plot(SNR(1:end),P_D_detection_estimation_GLRT(1:end),'b-o')
plot(SNR(1:end),P_D_detection_estimation_nonzeroCyclicAcf(1:end),'g-*')
plot(SNR(1:end),P_D_detection_estimation_cyclicCorr(1:end),'m-+')
legend('Proposed estimator','Based on [3]','Based on [2]')
xlabel('SNR')
ylabel('p_d')
ylim([0,1])
title(['Joint probability of detection and correctly estimating the cycle period (P=',num2str(P),', N=',num2str(N),', L=',num2str(L), ', M=',...
    num2str(M), ', D=',num2str(numD),')'])
