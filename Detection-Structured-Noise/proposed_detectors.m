function [GLRT, LMPIT] = proposed_detectors(z,L,P,noise_type)
%PROPOSED_DETECTORS 
%   This computes the GLRTs and LMPITs / LMPIT-inspired tests proposed in
%   A. Pries, D. Ramírez, P. J. Schreier "LMPIT-inspired Tests for Detecting
%   a Cyclostationary Signal in Noise with Spatio-Temporal Structure," (submitted)
%
%   Input:
%       z           -   received data
%       L           -   number of antennas
%       P           -   cycle period to be tested
%       noise_type  -   can be either (spatially) 'uncorrelated', 
%                       (temporally) 'white' or 'both'
%
%   Output:
%       GLRT        -   GLRT test statistic
%       LMPIT       -   the LMPIT if it exists, otherwise our
%                       LMPIT-inspired tests

if strcmp(noise_type,'uncorrelated') % spatially uncorrelated
    C_eig = eigCoherenceMatrix(z,1,L*P);
    
    GLRT = sum(log(C_eig));
    LMPIT = sum(C_eig.^2);
elseif strcmp(noise_type,'white') || strcmp(noise_type,'both')    
    N = size(z,1)/(L*P);
    M = size(z,2);
    
    % under H1, we obtain the eigenvalues of S by the singular values
    % of z - this is the same for both noise types
    % We use cellfun to vectorize this
    z_tilde = mat2cell(z,L*P*ones(1,N),M);
    eig_S1 = (cell2mat(cellfun(@svd,z_tilde,'UniformOutput',false))/sqrt((N*L*P))).^2;

    if strcmp(noise_type,'white')
        R_LP_blocks = cellfun(@mtimes,z_tilde,cellfun(@ctranspose,z_tilde,'UniformOutput',false),'UniformOutput',false)';
        R_LP_blocks = cellfun(@times, R_LP_blocks, repmat({1/(N*L*P)}, 1, N), 'UniformOutput', false);
        
        % estimate diagonal blocks of cov matrix and average them
        z_tilde = mat2cell(z/sqrt(N*L*P),L*ones(1,N*P),M);
        S_0_blocks1 = cellfun(@mtimes,z_tilde,cellfun(@ctranspose,z_tilde,'UniformOutput',false),'UniformOutput',false)';
        S_0_average = sum(cat(3,S_0_blocks1{:}),3)/(N*P);
        
        eig_S0 = eig(S_0_average);
        eig_S0(eig_S0 < 0) = 1e-10;
        eig_S0 = real(eig_S0); % evd could possibly be replaced by svd, then this does not happen
        
        GLRT = sum(log(eig_S1))-N*P*sum(log(eig_S0));

        % LMPIT does not exist, so return our LMPIT-inspired tests instead
        
        % This is the Frobenius norm of C_j
        norm_blocks = cellfun(@mtimes,R_LP_blocks,repmat({kron(eye(P),inv(S_0_average))},1,N),'UniformOutput',false);
        eig_new = cellfun(@eig, norm_blocks, 'UniformOutput',false);
        LMPIT.reg = sum(sum(abs(cell2mat(eig_new)).^2));

        % This is the Frobenius norm of C_av
        A = repmat({kron(eye(P),inv(sqrtm(S_0_average)))},1,N);
        C_hat_blocks = cellfun(@mtimes,cellfun(@mtimes,A,R_LP_blocks,'UniformOutput',false),A,'UniformOutput',false);
        C_hat_av = sum(cat(3,C_hat_blocks{:}),3)/N;
        LMPIT.av = norm(C_hat_av,'fro')^2;
    elseif strcmp(noise_type,'both')
        % estimate diagonal elements of cov matrix and average them
        S_0_blocks = sum(z.*conj(z),2)/(N*L*P);
        S_0_averaged = mean(reshape(S_0_blocks,L,N*P),2);

        GLRT = sum(log(eig_S1))-N*P*sum(log(S_0_averaged));
    end      
end