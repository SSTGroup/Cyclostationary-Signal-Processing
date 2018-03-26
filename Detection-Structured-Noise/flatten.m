function f = flatten(x)
%FLATTEN 
%   Reproduces np.flatten() which should be integrated in core Matlab
%
%   Input:
%       x   -   array
%
%   Output:
%       f   -   vector

f = reshape(x, [], 1);