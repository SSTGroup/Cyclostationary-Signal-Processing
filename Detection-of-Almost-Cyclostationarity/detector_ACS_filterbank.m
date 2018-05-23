function [statistic_GLRT, id_GLRT, frac_GLRT] = detector_ACS_filterbank(x,M, N_samples, P_max, numD, factor_matrix, N_len, bta)
% DETECTOR_ACS_FILTERBANK: Obtains the GLRT proposed in [1] for unknown 
% cycle period.
% To this end the signal is resampled with a filterbank as proposed in [2] 
% for every candidate integer part according to a given set of fractional 
% parts of the cycle period such that the resampled signal has an 
% integer-valued (or sufficiently small fractional part) cycle period. For 
% every candidate integer period we can apply the detector as proposed in 
% [1].
%
%           statistic_GLRT = detector_ACS_filterbank(RxSignal,M, N_samples, P_max, numD)
%           x = multivariate time series
%           M = number of snapshots (of length N_samples)
%           N_samples = number of samples per snapshot
%           P_max = maximum expected integer part of cycle period
%           numD = grid size of fractional parts of cycle period
%           factor_matrix = contains the prime factorization for each upsampling rate
%           N_len = length of kaiser window for the interpolation filter
%           bta = window parameter of kaiser window for interpolation filter
%
%           statistic_GLRT = detector statistic
%
% [1] D. Ramirez, P.J. Schreier, J. Via, I. Santamaria, and L.L. Scharf,
%     "Detection of multivariate cyclostationarity," IEEE Trans. Signal
%     Process., vol. 63, no. 20, p. 5395-5408, 2015
% [2] S. Horstmann, D. Ramirez, and P.J. Schreier, "Joint Detection of 
%     Almost-Cyclostationary Signals and Estimation of their Cycle Period" 
%     Submitted in IEEE Signal Processing Letters, May 2017

L = size(x,1); % Number of antennas
N = floor(N_samples*2/2.5); % Specify number of samples to consider after resampling (common for all P_int)

id_GLRT = zeros(P_max,1);
statistic_GLRT = zeros(P_max,1);
frac_GLRT = zeros(P_max,1);
%% Resampling stage -- filter bank implementation
% Upsample signal with L
epsilon = linspace(-0.5,0.5 - 1/numD,numD)';
p_vec = (2:P_max); p_vec = repmat(p_vec,numD,1); p_vec = reshape(p_vec,[],1);
p_frac_vec = p_vec + repmat(epsilon,P_max-1,1);
resam_rates_vec = p_vec./(p_frac_vec);

% Check if there are L/M=1 within the grid (no need for resampling)
id_1 = find(resam_rates_vec == 1);
[L_resam_orig,M_resam_orig] = rat(resam_rates_vec, 1e-12);

factor_matrix = sortrows(factor_matrix,'descend');

gamma_GLRT = zeros(numD*(P_max-1),1);


if ~isempty(id_1)
    p_vec_temp = ceil(id_1./numD) + 1;
    for ii = 1:length(id_1)
        x_reshape = reshape(x(:,1:floor(N/p_vec_temp(ii))*p_vec_temp(ii)*M),L,[],M);
        gamma_GLRT(id_1(ii)) = ...
            detector_cyclostationarity(x_reshape,p_vec_temp(ii),M);
    end
end

count_all = size(factor_matrix,1);
stage_matrix = factor_matrix;
signal_upsampled = cell(1,size(factor_matrix,2) + 1);
signal_upsampled(1) = {x.'};
ll = 1;
nn = ll;
temp_count = zeros(1,size(stage_matrix,2));
tic
while count_all >= 1
    if temp_count(ll) == 1
        while temp_count(ll) == 1
            index_upsam = find(L_resam_orig == prod(stage(1:ll)));
            p_vec_temp = ceil(index_upsam/numD)+1;
            for jj = 1:length(index_upsam)
                signal_resampled = downsample(signal_upsampled{ll+1},M_resam_orig(index_upsam(jj)));
                signal_resampled = reshape(signal_resampled(1:floor(N/p_vec_temp(jj))*p_vec_temp(jj)*M,:).',L,[],M);
                gamma_GLRT(index_upsam(jj)) = detector_cyclostationarity(...
                    signal_resampled(:,1:floor(N/p_vec_temp(jj))*p_vec_temp(jj),:),p_vec_temp(jj),M);
                clear signal_resampled
            end
            temp_count = temp_count - 1;
            stage_matrix = stage_matrix(2:end,:);
            count_all = count_all - 1;
            idx = find(temp_count(1:ll-1) == 0);
            if ~isempty(idx)
                ll = idx(1);
                nn = ll;
            else
                nn = ll;
                ll = ll - 1;
                if ll == 0
                    break
                end
            end
        end
        if ll == 0
            ll = 1;
        end
    else
        stage = stage_matrix(1,(stage_matrix(1,:) ~= 0));
        for mm = nn : length(stage)
            L_temp = stage(mm);
            fc = 1/2/L_temp;
            window_length_kaiser = 2*N_len*L_temp + 1;
            h_filt = firls(window_length_kaiser-1, [0 2*fc 2*fc 1], [1 1 0 0]).*kaiser(window_length_kaiser,bta)' ;
            h_filt = L_temp*h_filt/sum(h_filt);
            
            firinterp = dsp.FIRInterpolator(L_temp,h_filt);
            signal_upsampled(mm+1) = {firinterp(signal_upsampled{mm})};
            
            kk = 0;
            temp_count(mm) = 0;
            while temp_count(mm) == 0
                if sum(stage_matrix(kk+1,1:mm) == stage(1:mm)) == mm
                    kk = kk + 1;
                    if kk + 1 > size(stage_matrix,1)
                        temp_count(mm) = kk;
                    end
                else
                    temp_count(mm) = kk;
                end
            end
        end
        ll = mm;
    end
end

for p = 1:P_max-1
    [statistic_GLRT(p+1),id_GLRT(p)] = min(gamma_GLRT(numD*(p-1)+1:numD*p)); % Find maximum likelihood estimate of fractional part by minimizing the GLRT statistic
    frac_GLRT(p+1) = epsilon(id_GLRT(p));
end
end