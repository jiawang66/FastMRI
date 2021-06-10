function krecov = pocs_2d(kdata, pff, dim, delta, maxs)
% function krecov = pocs(kdata)
%   Projection Onto Convex Sets (POCS), a iterative methods for partial
%   Fourier reconstruction of 2-D magnetive resonance imaging (MRI).
%
% Input -
% kdata: 2-D k-space
% pff: partial Fourier fraction, 0.5 <= pff <= 1
% dim: partial sample along one dimension, 1 <= dim <= N
% delta: iteration termination condition, default delta = 1
% maxs: maximum iteration steps, default maxs = 100
%
% Output -
% krecov: recovered k-space

%% Parameters check
if nargin < 3 || pff < 0.5 || pff > 1
    error("Missing parameters or wrong 'pff'. (0.5 <= pff <= 1) needed");
end

N = numel(size(kdata));
if N > 2
    error("2-D data expected.")
end

if dim < 1 || dim > N
    error("Wrong 'dim'. (1 <= dim <= N) needed.");
end

if nargin < 4
    delta = 1;    % Iteration termination condition
    maxs = 1e2;     % maximum iteration steps
elseif nargin < 5
    maxs = 1e2;
end

%% Initialization
kdata = autopermute(kdata, dim);
up = round(size(kdata, 1) * pff) - 1;   % cutoff_up
down = size(kdata, 1) - up + 1;         % cutoff_down, about k-space's center symmetry

% cosine cone window
hw = 5;     % half width
f = @(x,c,w) (cos(pi/(4*w) * x + pi/4 - pi*c/(4*w))).^2;

win = zeros(size(kdata,1),1);
win(1:up) = 1;
if up + hw <= size(kdata,1)
    win((up-hw):(up+hw)) = f((up-hw):(up+hw), up, hw);
end
win = repmat(win, 1, size(kdata,2));

% prepare matrix
ksym = zeros(size(kdata));
krecov = zeros(size(kdata));
ktemp = zeros(size(kdata));

ksym(down:up, :) = kdata(down:up, :);   % symmetric data
krecov(1:up, :) = kdata(1:up, :);       % constraint data

%% IDFT of the symmetric data
msym = idft2(ksym);
phase = msym ./ abs(msym);

%% Iterative reconstruction
count = 0;
bias = matnorm(krecov(1:up,:) - ktemp(1:up,:), 3);
disp(['relative error : ', num2str(bias)])

while count < maxs && bias > delta
    count = count + 1;
    mrecov = abs(idft2(krecov)) .* phase;
    ktemp = dft2(mrecov);
    krecov = win .* krecov + (1 - win) .* ktemp;
%     krecov(up+1:end, :) = ktemp(up+1:end, :);
    bias = matnorm(krecov(1:up,:) - ktemp(1:up,:), 3);
    disp(['relative error : ', num2str(bias)])
end

krecov = autopermute(krecov, dim);

end