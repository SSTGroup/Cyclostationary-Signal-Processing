function [gamma_GLRT,gamma_LMPIT_LSR,varargout] = passive_detector_cyclostationarity(u_s,u_r,P,M)
% PASSIVE_DETECTOR_CYCLOSTATIONARITY: Obtains the GLRT and an LMPIT-
% inspired test for the two-channel passive detection problem exploiting 
% cyclostationarity.
%
%           [gamma_GLRT,gamma_LMPIT_LSR] = passive_detector_cyclostationarity(u,P,M)
%           u_s = multivariate time series observed at surveillance channel
%           u_r = multivariate time series observed at reference channel
%           P = cycle period
%           M = number of snapshots (of length L*P*N)
%
%           Output: 
%           gamma_GLRT = GLRT statistic
%           gamma_LMPIT_LSR = LMPIT-inspired test statistic
%
%           Optional output arguments:
%           gamma_LMPIT_LS = LMPIT statistic for detecting a CS signal at 
%                           the  surveillance channel only
%           gamma_GLRT_S = GLRT statistic for detecting a CS signal at the 
%                           the surveillance channel only
%
% Assume an equal number of channels at SC and RC
[L,NMP] = size(u_s);

% Warning if we're getting singular matrices!
if M < 2*L*P
    warning('M < 2LP! Matrices will be singular!')
end
N = NMP/(P*M);

% Divide long observation vector into M windows and treat them as if they
% were i.i.d. observations. 
w_s = reshape(u_s,L*N*P,M);
w_r = reshape(u_r,L*N*P,M);

% Transform vector-valued sequences to the Fourier-Domain 
x = reshape(ifft(reshape(w_s,L,N*P,M),[],2),L*N*P,M);
y = reshape(ifft(reshape(w_r,L,N*P,M),[],2),L*N*P,M);

% Obtain the sample (cross-) covariance matrices
XX = x*x'/M;
YY = y*y'/M;
XY = x*y'/M;

% Obtain commutation matrix
L_NP_N = zeros(N*P);
for idx = 0:P-1
    for jdx = 0:N-1
        L_NP_N(jdx*P+idx+1,idx*N+jdx+1) = 1;
    end
end
L_kron = kron(L_NP_N,eye(L));

% Reorder elements of the sample (cross-) covariance matrices to obtain the
% blocks of the sample covariance matrix of z (as given in the manuscript).
Qs = L_kron*XX*L_kron';
Qr = L_kron*YY*L_kron';
Qsr = L_kron*XY*L_kron';

%% Obtain the test statistics
gamma_LMPIT_LSR = 0;
gamma_LMPIT_LS = 0;
gamma_GLRT_SR = 0;
gamma_GLRT_S = 0;

for k = 1:N
    % Get the LP x LP sized blocks of the main diagonal of Qs, Qr, and Qsr,
    % respectively
    Qs_temp = Qs((k-1)*L*P+1:k*L*P,(k-1)*L*P+1:k*L*P);
    Qr_temp = Qr((k-1)*L*P+1:k*L*P,(k-1)*L*P+1:k*L*P);
    Qsr_temp = Qsr((k-1)*L*P+1:k*L*P,(k-1)*L*P+1:k*L*P);
    
    % Get the L x L sized blocks of the main diagonal of Qs
    diag_L_Qs_temp = zeros(L*P);
    for ii = 1:P
        diag_L_Qs_temp((ii-1)*L+1:ii*L,(ii-1)*L+1:ii*L) = Qs_temp((ii-1)*L+1:ii*L,(ii-1)*L+1:ii*L);
    end
    
    % Obtain L_S as given in Theorem 2 in the manuscript
    C11_temp = sqrtm(diag_L_Qs_temp)\Qs_temp/sqrtm(diag_L_Qs_temp);
    svd_C11 = real(svd(C11_temp));
    gamma_LMPIT_LS = gamma_LMPIT_LS + sum(svd_C11.^2);
    
    % Obtain the LMPIT-inspired test L_SR as given in Theorem 2 in the 
    % manuscript
    C12_temp = sqrtm(diag_L_Qs_temp)\Qsr_temp/sqrtm(Qr_temp);
    gamma_LMPIT_LSR = gamma_LMPIT_LSR + sum(real(svd(C12_temp)).^2);
    
    % Obtain the GLRT statistic
    gamma_GLRT_S = gamma_GLRT_S + sum(log(svd_C11));
    
    gamma_GLRT_SR = gamma_GLRT_SR + sum(log(real(eig(Qr_temp - Qsr_temp'/Qs_temp*Qsr_temp,Qr_temp)))); 
end

gamma_GLRT = gamma_GLRT_S + gamma_GLRT_SR;

varargout{1} = gamma_LMPIT_LS;
varargout{2} = gamma_GLRT_S;
