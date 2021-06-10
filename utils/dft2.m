function y = dft2(x)
% function f = dft2(x)
%   2-Dimension DFT for 2-D matrix through twice 1-D DFT

%% Parameters check
if numel(size(x)) > 2
    error("2-D data Expected");
end

%% Compute
y = dftmatrix(size(x,1), 1) * double(x);
y = dftmatrix(size(x,2), 1) * y';
y = y';
end