function y = dftn(x, n, dims)
% function y = dftn(x, n, dims)
%   This MATLAB function returns the multidimensional DFT of N-D array in 
%   arbitrary dimensions (given by parameter 'dims') when n = -1. When n >
%   0, n times DFT will do successively in first n dimensions.

%% Parameters check
if n == -1 && nargin < 3
    dims = 1 : numel(size(x));
elseif n > 0
    if numel(size(x)) < n
        error("Wrong dimensions");
    end
    dims = 1 : n;
end

%% compute
x = double(x);

for k = dims
    temp = autopermute(x, k);
    shape = size(temp);
    temp = reshape(temp, shape(1), []);
    x = dftmatrix(shape(1), 1) * temp;
    x = autopermute(reshape(x, shape), k);
end

y = x;

end