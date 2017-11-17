function gamma_GLRT = detector_cyclostationarity(u,P,M)
% DETECTOR_CYCLOSTATIONARITY: Obtains the GLRT for the test cyclostationarity
% with period P vs wide-sense stationarity
%
%
%           [gamma_GLRT,gamma_LMPIT] = detector_cyclostationarity(u,P,M)
%           u = multivariate time series
%           P = cycle-period
%           M = number of snapshots (of length P*N) 
%       
%           gamma_GLRT = GLRT statistic

% David Ramirez
% Signal and System Theory group
% University of Paderborn (Germany)
% 2014

[L,NMP] = size(u);
N = NMP/(P*M);

%% Obtain the covariance matrix

X = reshape(u,L*N*P,M);
Rxx = X*X'/M;

%% Pre- and post-multiply by the Fourier matrix

dummy = reshape(ifft(reshape(Rxx,L,N*P,L*N*P),[],2),L*N*P,L*N*P);
Sxx = reshape(ifft(reshape(dummy',L,N*P,L*N*P),[],2),L*N*P,L*N*P)*P*N;

%% Commutation matrix

K = sparse([],[],[],N*P,N*P,N*P);

K(N*P,N*P) = 1;
index_2 = 0:N*P-2;
index_1 = mod(index_2*N, N*P-1);
K(index_1*N*P+index_2+1) = 1;

K_bar = kron(K,speye(L));

Sxx = K_bar*Sxx*K_bar';

%% Obtain the eigenvalues of the Coherence matrix and the statistics

gamma_GLRT = 0;

mask = logical(kron(speye(P),ones(L,L)));

for k = 1:N
    S_k = Sxx((k-1)*L*P+1:k*L*P,(k-1)*L*P+1:k*L*P);
    D_k = S_k;D_k(~mask) = 0;
    eigenvalues = abs(real(eig(S_k,D_k)));  % In case there are complex eigenvalues with negative real part due to numerical problems
    gamma_GLRT = gamma_GLRT + sum(log(eigenvalues));
end