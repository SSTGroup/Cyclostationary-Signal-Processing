function n = structured_noise(L,N_samples,spatial_correlation,temp_white)
%STRUCTURED_NOISE 
%   Convenience function to sample noise with certain spatial or temporal
%   structure/correlation.
%
%   Input:
%       L           -   number of antennas
%       N_samples   -   number of samples
%       spatial_correlation   -   matrix to add spatial correlation
%       temp_white  -   whether to sample temporally white noise (bool)
%
%   Output:
%       n           -   resulting noise signal

n = randn(L,N_samples)+ 1j*randn(L,N_samples);
if ~temp_white
    % filtering is only needed if noise is not white
    n = filter(1/19*ones(1,19),1,n,[],2); % moving average
end
% Add spatial correlation
n = spatial_correlation*n;