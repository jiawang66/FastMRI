function y = normabs(x, method)
% function y = normabs(x)
%   compute the abs after normalization

%% Preparation
if nargin < 2
    method = 'range';
end

%% compute
shape = size(x);
y = normalize(reshape(abs(x), numel(x), []), method);
y = reshape(y, shape);
end