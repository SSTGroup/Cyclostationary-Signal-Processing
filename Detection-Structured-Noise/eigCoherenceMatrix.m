function eig_C = eigCoherenceMatrix(z,B0,B1)
%EIGCOHERENCEMATRIX Obtain eigenvalues of coherence matrix    
%   The block-diagonal coherence matrix is defined by 
%       diag_B0 (S)^{-0.5} [diag_B1 (S)] diag_B0 (S)^{-0.5}
%
%   Input:
%       S   -   covariance matrix
%       B0  -   outer block size
%       B1  -   inner block size (B1 > B0)
%
%   Output:
%       C   -   array with eigenvalues

B_ratio = B1/B0;
mask_B0 = kron(eye(B_ratio),ones(B0));

eig_C = zeros(size(z,1)/B1,1);
for b=1:size(z,1)/B1
    % estimate part of covariance matrix with big block size B1
    D1 = (z((b-1)*B1+1:b*B1,:)*z((b-1)*B1+1:b*B1,:)')/(size(z,1));
    
    % estimate part of covariance matrix with small block size B0
    % by reusing parts of D1
    D0 = zeros(B1);
    D0(mask_B0==1) = D1(mask_B0==1);
    
    % obtain eigenvalues of coherence matrix by the generalized
    % evd
    eig_C((b-1)*B1+1:b*B1) = eig(D1,D0);
end

%% integrated self-checks / numerical stability
% all eigenvalues should be real, so print a warning if that's not the case
if max(abs(imag(eig_C))) > 1e-5
    % give warning, if imaginary part is not small
    warning('Imaginary part')
end
eig_C = real(eig_C);

% Make sure that small negative numbers are set to a small
% positive number. This avoids numerical issues if we compute log(eig)
if min(eig_C) < -1e-5
    % give warning, if not small negative number
    warning('Imaginary part')
end
eig_C(eig_C < 0) = 1e-10;

end