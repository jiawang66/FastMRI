function y = autopermute(x, dim)
% function y = autopermute(x, dim)
%   Auto exchange the dimension 1 and 'dim' of N-dimension array x

%% Parameters check
if nargin < 2
    error('Parameter "dim" needed.');
end

N = numel(size(x));     % number of dimensions

if dim < 1 || dim > N || round(dim) ~= dim
    error('Expected 1 < dim < N and dim should be an integer.');
end

%% exchange
dims = 1 : N;
dims(dim) = 1;
dims(1) = dim;

y = permute(x, dims);

end