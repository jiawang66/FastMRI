function y = myfftshift(x ,dim)
% Function overloading of inner function fftshift(), myfftshift() do
% "half-spaces" swaping along any dimension specified by parameter vector
% dim.

%% prepare
if nargin < 2
    dim = 1 : numel(size(x));
end

%% swap
y = x;
for k = 1 : numel(dim)
    y = fftshift(y, dim(k));
end

end