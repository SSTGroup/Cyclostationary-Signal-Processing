function [gamma_GLRT] = passive_detector_corr_subspace(u_s,u_r,rho)
% PASSIVE_DETECTOR_CORR_SUBSPACE: Obtains the GLRT for the detection of 
% correlated subspaces in two MIMO channels for noise model 4 as  proposed 
% in [1] detection problem.
%
%           [gamma_GLRT] = passive_detector_corr_subspace(u_s,u_r,rho)
%           u_s = multivariate time series observed at surveillance channel
%           u_r = multivariate time series observed at reference channel
%           rho = number of antennas at the IO
%
%           Output: 
%           gamma_GLRT = GLRT statistic
%
[~,NMP] = size(u_s);

% Obtain sample (cross-) covariance matrices and the cohrence matrix and 
% its SVD
S_ss = u_s*u_s'/NMP;
S_sr = u_s*u_r'/NMP;
S_rr = u_r*u_r'/NMP;

C = sqrtm(inv(S_ss))*S_sr*sqrtm(inv(S_rr));
k = svd(C);
k = sort(k,'descend');

% Obtain the test statistic
gamma_prod = zeros(rho,1);
for ii = 1:rho
    gamma_prod(ii) = 1/(1-k(ii)^2);
end
gamma_GLRT = prod(gamma_prod);
end
