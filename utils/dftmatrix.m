function DM = dftmatrix(N, f)
% function DM = dftmatrix(N, f)
%   Generate NxN DFT (Discrete Fourier Transform) matrix for 1-D FT
%
% Input -
% N: Size of the DFT matrix
% f: 1 -- DFT, 2 -- IDFT, default f = 1
%
% Output -
% DM: NxN DFT matrix

%% Parameter check
if N ~= round(N) || N < 1
    error('N needed to be a positive integer (N > 1)')
end

if nargin < 2
    f = 1;
elseif f ~= 1 && f ~= 2
    error('Inputed parameter "f" should be 1 or 2.')
end

if f == 2
    f = -1;
end

%% Generating
w = exp(-1i * 2 * pi / N * f);
DM = ones(N);

alpha = (1 : N-1)';
subM = alpha * alpha';

DM(2:N, 2:N) = w .^ subM;

if f == -1
    DM = DM / N;
end
end