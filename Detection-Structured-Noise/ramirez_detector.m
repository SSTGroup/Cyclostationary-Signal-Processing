function [GLRT,LMPIT] = ramirez_detector(z,L,P)
%RAMIREZ_DETECTOR GLRT and LMPIT detection of cyclostationarity
% 
% Based on the journal paper
%   D. Ramirez, P. J. Schreier, J. Vía, I. Santamaria, and L. L. Scharf, 
%   "Detection of multivariate cyclostationarity," IEEE Trans. Signal 
%   Process., vol. 63, no. 20, pp. 5395?5408, 2015.
%
%   Input:
%       z   -   columns are individual realizations of z
%       L   -   number of receiver antennas
%       P   -   samples per cycle period
%
%   Output:
%       GLRT    -  GLRT test statistic
%       LMPIT   -  LMPIT test statistic
%

% GLRT
C_eig = eigCoherenceMatrix(z,L,L*P);
GLRT = sum(log(C_eig));

% LMPIT
LMPIT = sum(C_eig.^2);

end