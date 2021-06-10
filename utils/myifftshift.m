function y = myifftshift(x ,dim)
% Function overloading of inner function ifftshift(), myifftshift() do
% "half-spaces" swaping along any dimension specified by parameter vector
% dim.

%% prepare
if nargin < 2
    dim = 1 : numel(size(x));
end

%% swap
y = x;
for k = 1 : numel(dim)
    y = ifftshift(y, dim(k));
end

end