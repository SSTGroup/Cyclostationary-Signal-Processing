function z = getAllZ(N, P, M, x)
%GETALLZ 
%   Convenience function to compute a set of realizations from a single
%   realization.
%
%   Input:
%       N, P, N   -   signal parameters length, cycle period, and
%       realizations
%       x         -   matrix-valued signal (L x MNP)
%
%   Output:
%       f   -   vector

L = size(x,1);
assert(N*P*M == size(x,2), 'invalid dimensions')

z = zeros(numel(x)/M, M);
% Chop the signal into M pieces of length LNP
for m = 1:M
    m_ind = (m-1)*N*P + [1:N*P];
    z(:, m) = getZ(L, N, P, x(:,m_ind));
end