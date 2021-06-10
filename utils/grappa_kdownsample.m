function [dkdata, data] = grappa_kdownsample(kdata, rate, dim)
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
data = reshape(data, shape(1), []);

% --- mod(length_of_the_1st_dimension, rate) should be 1, then the first
% and last line of downsample data is not empty
if mod(shape(1), rate) ~= 1
    if mod(shape(1), rate) == 0
        new_line = shape(1) - rate + 1;
    else
        new_line = shape(1) - mod(shape(1), rate) + 1;
    end

    new_start = floor((shape(1) - new_line) / 2) + 1;
    new_indices = new_start : (new_start + new_line - 1);
    data = data(new_indices, :);
    shape(1) = new_line;
end

% --- downsample
dkdata = zeros(size(data));
dkdata(1:rate:end, :) = data(1:rate:end, :);
dkdata = reshape(dkdata, shape);
dkdata = autopermute(dkdata, dim);

data = reshape(data, shape);
data = autopermute(data, dim);

end