function z = getZ(L, N, P, x)
%GETZ Computes the transformation x -> z
%   This computes an efficient version of Eq. (8) in our paper. To avoid
%   the expensive matrix multiplications, we compute the first part by an
%   FFT (actually an IFFT). The second product does not need any FLOPs as
%   it is a simple re-ordering of elements. Everything is vectorized.
%
%   Input:
%       L           -   number of channels
%       h_profile   -   power delay profile (vector)
%
%   Output:
%       h_channel   -   channel taps

% To avoid computing the indexing every time in a Monte Carlo simulation, 
% we cache the result of this function. This creates issues when this is 
% called with different signal parameters and requires calling 
% `clear functions` to force re-computation.
persistent j_sorted;

if isempty(j_sorted)
    % If the indexing is not cached yet compute it
    [i,j] = find(kron(commutationMatrix(P,N),eye(L)));
    [~,i_sorted] = sort(i);
    j_sorted = j(i_sorted);
end

if numel(j_sorted)/L ~= size(x,2)
    % throw an error if an incorrect indexing is still cached
    error('Unexpected error, %d, %d', numel(j_sorted), size(x,2))
end

% Actual computation
x_ifft = N*P*ifft(x,[],2);
z = x_ifft(j_sorted);