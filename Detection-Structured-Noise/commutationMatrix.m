function C = commutationMatrix(size_A,varargin)
%COMMUTATIONMATRIX Obtain a commutation matrix of specific size
%   commutation matrix C, such that for any matrix A with size(A) = size_A
%   A(:) = C*A'(:)
%
%   Input:
%       size_A   -   two-element vector
%
%   Output:
%       C        -   commutation matrix
%

if nargin == 2
    size_A = [size_A varargin{:}];
end

C_helper = reshape(1:prod(size_A),size_A)';
C = zeros(prod(size_A));
C(sub2ind(size(C),(1:size(C,1))',C_helper(:))) = 1;
C = C';

% test with this code-snippet
    % A = randn(size_C);
    % AT = A';
    % norm(A(:)-C*AT(:))
end

