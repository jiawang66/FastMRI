function dkdata = kdownsample(kdata, r, dim)
% function d = kdownsample(kdata, r)
%   downsample by factor 'r' in dimension 'dim'
%   Assume that k-space center is in the center of the matrix kdata.

%% Parameters check
if dim < 1 || dim > numel(size(kdata)) || round(dim) ~= dim
    error('Expected parameter dim to be an integer.')
end

%% Downsample
% --- permute first dimension with dimension dim
data = autopermute(kdata, dim);
shape = size(data);
data = reshape(data, shape(1),[]);

% --- downsample
dkdata = zeros(size(data));
dkdata(1:r:end, :) = data(1:r:end, :);
dkdata = reshape(dkdata, shape);
dkdata = autopermute(dkdata, dim);

end