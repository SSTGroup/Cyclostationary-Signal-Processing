function [gamma_GLRT,gamma_LMPIT] = detectors_Gardner(x,x_shifted,tau_0)
% DETECTORS_GARDNER: Obtains the GLRT and LMPIT for a fixed lag and cycle
% frequency
%
%
%           [gamma_GLRT,gamma_LMPIT,gamma_Cabric] = detectors_fixed_lag(x,x_shifted,tau_0)
%           x = vector-valued time series
%           x_shifted = frequency-shifted version of x
%           tau_0 = time lag used for the covariance matrices
%       
%           gamma_GLRT = GLRT statistic
%           gamma_LMPIT = LMPIT statistic

% David Ramirez
% Signal and System Theory group
% University of Paderborn (Germany)
% 2014

if ndims(x) == 3
    L = size(x,1);
    x = reshape(x,L,[]);
    x_shifted = reshape(x_shifted,L,[]);
end
N = size(x,2);

%% Covariance matrices

Rxx = x*x'/N;

%% Cyclic covariance matrices

Rxx_alpha = x_shifted(:,1+tau_0:end)*x(:,1:end-tau_0)'/(N-tau_0);

%% Coherence matrix and svd

[V,D] = eig(Rxx);
Rxx_12 = V*diag(diag(D).^(-1/2))*V';

C = Rxx_12*Rxx_alpha*Rxx_12;

mu = svd(C);

%% GLRT

gamma_GLRT = -sum(log(1-mu.^2));

%% LMPIT

gamma_LMPIT = -sum(mu.^2);
