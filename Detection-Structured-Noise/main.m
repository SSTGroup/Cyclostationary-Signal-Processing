% This is a MATLAB implementation supporting the paper
%
% "LMPIT-inspired Tests for Detecting a Cyclostationary Signal in Noise with
% Spatio-Temporal Structure" by Aaron Pries, David Ramirez and
% Peter J. Schreier, IEEE Transactions on Wireless Communications,
% vol. 17, no. 9, pp. 6321-6334, Sept. 2018
%
% ## ----------------------------------------------------------------------------
% ##
% ##   File: main.m
% ##   Copyright (c) <2018> Signal and System Theory Group,
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
% ##   Author: Aaron Pries
% ##   Date: 09.03.2018
% ##   Version: v1
% ##
% ##   At last compiled with Matlab 2016b
% ##   Required Toolboxes:
% ##        Communications System Toolbox
% ##        Signal Processing Toolbox
% ##        Statistics Toolbox
% ## ----------------------------------------------------------------------------


clear variables; close all force; clc
clear functions

%% Description
% This script runs repeated experiments at different SNRs to compare the
% performance of different detectors of cyclostationarity. It recreates
% Figure 8 in our paper.

%% User Parameters
N_iter = 1000; % number of Monte Carlo experiments
N_SNR = 10; % how many SNRs we want to generate
SNR_list = linspace(-20,-5,N_SNR);

NM = 2^10;  % total observation time
M = 64; % number of 'realizations'
N = NM/M; % observation length per realization
L = 2; % number of antennas

% OFDM part (using the Communications Toolbox)
ofdmMod = comm.OFDMModulator('FFTLength',16, 'NumTransmitAntennas', L , ...
    'NumSymbols', NM, 'CyclicPrefixLength', 4, 'NumGuardBandCarriers', [0; 0], ...
    'InsertDCNull', false);

% Calculate the cycle period
P = ofdmMod.CyclicPrefixLength + ofdmMod.FFTLength;
ofdmModDim = info(ofdmMod);
numData = ofdmModDim.DataInputSize(1);

% Create a QPSK modulator.
qpskMod = comm.QPSKModulator;

T = 30; % channel taps
h_profile = sqrt(exp(-0.1*(0:T-1)));  % channel delay profile

% finally, check validity of parameters
if M < L*P
    error('Not enough realizations')
end

%% Monte Carlo Experiments
% allocate everything to be faster
LMPIT_0 = zeros(N_SNR,N_iter);LMPIT_1 = zeros(N_SNR,N_iter);
LMPIT_av_0 = zeros(N_SNR,N_iter);LMPIT_av_1 = zeros(N_SNR,N_iter);
GLRT_0 = zeros(N_SNR,N_iter);GLRT_1 = zeros(N_SNR,N_iter);
LMPIT_WSS_0 = zeros(N_SNR,N_iter);LMPIT_WSS_1 = zeros(N_SNR,N_iter);
ref_detector0 = zeros(N_SNR,N_iter);ref_detector1 = zeros(N_SNR,N_iter);
urriza_H0 = zeros(N_SNR,N_iter);urriza_H1 = zeros(N_SNR,N_iter);

for iter = 1:N_iter
    % Every Monte Carlo simulation uses a different channel and spatial
    % correlation
    h_channel = getLChannels(L, h_profile);
    sigma_w = randn(L)+1j*randn(L);

    % Allocate some variables
    z_s = zeros(L*N*P,M);
    z_w = zeros(L*N*P,M);
    LMPIT_av_1_temp = zeros(1,N_SNR);
    LMPIT_av_0_temp = zeros(1,N_SNR);
    GLRT_1_temp = zeros(1,N_SNR);
    GLRT_0_temp = zeros(1,N_SNR);
    LMPIT_WSS_1_temp = zeros(1,N_SNR);
    LMPIT_WSS_0_temp = zeros(1,N_SNR);
    ref_detector0_temp = zeros(1,N_SNR);
    ref_detector1_temp = zeros(1,N_SNR);
    urriza_H0_temp = zeros(1,N_SNR);
    urriza_H1_temp = zeros(1,N_SNR);

    for jj = 1:numel(SNR_list)
        SNR = SNR_list(jj); % SNR of the current experiment

        % Sample signal and noise - the same for all detectors
        w_long = structured_noise(L,N*P*M,sigma_w,true);
        qpsk_symbols = reshape(qpskMod(flatten(randi([0 3],N*numData*M, L))), numData, N*M, L);
        dataOFDM = ofdmMod(qpsk_symbols);
        s_long = multipleOFDM_signals(N*M, P, dataOFDM, h_channel);
        s_long = s_long*sqrt(10^(SNR/10))*norm(w_long,'fro')/norm(s_long,'fro');

        % For our detectors we assume independent realizations - here we
        % create them by chopping the long signal into fake realizations
        % for a fair comparison.
        z0 = getAllZ(N, P, M, w_long);
        z1 = getAllZ(N, P, M, s_long+w_long);

        % our GLRT and LMPIT-inspired detectors
        [GLRT_0_temp(jj), LMPIT_0_temp(jj)] = proposed_detectors(z0,L,P,'white');
        [GLRT_1_temp(jj), LMPIT_1_temp(jj)] = proposed_detectors(z1,L,P,'white');

        % LMPIT from Ramirez et al, 201
        [~,LMPIT_WSS_0_temp(jj)] = ramirez_detector(z0,L,P);
        [~,LMPIT_WSS_1_temp(jj)] = ramirez_detector(z1,L,P);

        % For the Lunden test the tests statistic is aggregated for
        % all antennas
        for l = 1:L
            % We use the first cycle frequency 1/P
            [~,temp0] = cycleTest(w_long(l,:),0,[1, -1] * ofdmMod.FFTLength,N*M+1);
            ref_detector0_temp(jj) = ref_detector0_temp(jj) + temp0;

            [~,temp1] = cycleTest(s_long(l,:)+w_long(l,:),0,[1, -1] * ofdmMod.FFTLength,N*M+1);
            ref_detector1_temp(jj) = ref_detector1_temp(jj) + temp1;
        end

        % Evaluate the Urriza detector at the first cycle frequency 1/P
        urriza_H0_temp(jj) = urriza_detector(w_long,-ofdmMod.FFTLength,1/P,1);
        urriza_H1_temp(jj) = urriza_detector(s_long+w_long,-ofdmMod.FFTLength,1/P,1);
    end
    LMPIT_1(:,iter) = cell2mat({LMPIT_1_temp.reg});
    LMPIT_0(:,iter) = cell2mat({LMPIT_0_temp.reg});
    LMPIT_av_1(:,iter) = cell2mat({LMPIT_1_temp.av});
    LMPIT_av_0(:,iter) = cell2mat({LMPIT_0_temp.av});
    GLRT_1(:,iter) = GLRT_1_temp;
    GLRT_0(:,iter) = GLRT_0_temp;
    LMPIT_WSS_1(:,iter) = LMPIT_WSS_1_temp;
    LMPIT_WSS_0(:,iter) = LMPIT_WSS_0_temp;
    ref_detector0(:,iter) = ref_detector0_temp;
    ref_detector1(:,iter) = ref_detector1_temp;
    urriza_H0(:,iter) = urriza_H0_temp;
    urriza_H1(:,iter) = urriza_H1_temp;
end

%% Evaluation: obtain P_MD and plot

P_FA = 0.01; % Probability of false alarm
for kk = 1:numel(SNR_list)
    GLRT_PD(kk) = sum(GLRT_1(kk,:) < prctile(GLRT_0(kk,:),P_FA*100))/N_iter;
    LMPIT_PD(kk) = sum(LMPIT_1(kk,:) > prctile(LMPIT_0(kk,:),100-P_FA*100))/N_iter;
    LMPIT_av_PD(kk) = sum(LMPIT_av_1(kk,:) > prctile(LMPIT_av_0(kk,:),100-P_FA*100))/N_iter;
    LMPIT_WSS_PD(kk) = sum(LMPIT_WSS_1(kk,:) > prctile(LMPIT_WSS_0(kk,:),100-P_FA*100))/N_iter;
    urriza_PD(kk) = sum(urriza_H1(kk,:) > prctile(urriza_H0(kk,:),100-P_FA*100))/N_iter;
    lunden_PD(kk) = sum(ref_detector1(kk,:) > prctile(ref_detector0(kk,:),100-P_FA*100))/N_iter;
end

figure,hold all
pm4 = plot(SNR_list,LMPIT_WSS_PD,'+-','Color',[0 .6 0],'DisplayName','1','LineWidth',1.25,'MarkerSize',7);
pm1 = plot(SNR_list,GLRT_PD,'b*-','DisplayName','2','LineWidth',1.25,'MarkerSize',7);
pm2 = plot(SNR_list,LMPIT_PD,'rd-','DisplayName','3','LineWidth',1.25,'MarkerSize',4);
pm3 = plot(SNR_list,LMPIT_av_PD,'kx-','DisplayName','4','LineWidth',1.25,'MarkerSize',7);
pm5 = plot(SNR_list,lunden_PD,'>-','Color',[0 .7 .7],'DisplayName','Lunden et al.','LineWidth',1.25,'MarkerSize',7);
pm6 = plot(SNR_list,urriza_PD,'m<-','DisplayName','Urriza et al.','LineWidth',1.25,'MarkerSize',7);

ylabel('$P_{\mathrm{D}}$','Interpreter','Latex')
xlabel('$\mathrm{SNR}$','Interpreter','Latex')
legend([pm4,pm1,pm2,pm3,pm5,pm6],'Location','NorthWest')
