function recover_kspace = grappa_2d(downsample_kspace, rate, dim, kernel_size)

% function recover_kspace = grappa_2d(downsample_kspace, rate, dim, kernel_size)
%
%   GeneRalized Autocalibrating Partially Parallel Acquisitions (GRAPPA)
%   reconstruction for 2-D MR image, which is in 3-D matrix of data format
%   with the dimensions [pe, ro, channel] or [ro, pe, channel].
%
% Input --
% downsample_kspace: Cartisian acquired k-space data with zero-fill in the
% unacquired data line.
% rate: acceralation rate, must be less than the number of channels
% dim: scalar number, acceralation dimension
% kernel_size: two-element matrix, the size of GRAPPA kernel
%
% Output --
% recover_kspace: k-space recovered by GRAPPA method, default [2, 3]

%% Parameters preparation
assert(numel(size(downsample_kspace)) == 3, 'The inputed data should be 3-D [pe, ro, channel]')
assert(size(downsample_kspace,3) > rate, 'Acceralation rate should be less than number of channels')
assert(numel(size(downsample_kspace)) > dim, 'Wrong dimension. Prefer dim = 1 or dim = 2 for 2-D MRI')

% --- set default kernel size
if nargin < 4
    kernel_size = [2, 3];
end

assert(ismatrix(kernel_size) && numel(kernel_size)==2, 'Kernel_size should be a two-elememt matrix')
assert(kernel_size(1) > 1, 'First parameter of kernel size should larger than 1.')

%% Initialization
downsample_kspace = autopermute(downsample_kspace, dim);
[num_phase_orig, num_freq, num_coil] = size(downsample_kspace);

% --- Expect that the first and last line along phase dimension were
% sampled. If not, zero padding at the end of the phase dimension to ensure
% that mod(num_phase, rate) = 0
if mod(num_phase_orig, rate) ~= 1
    if mod(num_phase_orig, rate) ~= 0
        num_phase = num_phase_orig + rate - mod(num_phase_orig, rate);
        tmp = zeros([num_phase, num_freq, num_coil]);
        tmp(1:num_phase_orig, :, :) = downsample_kspace;
        downsample_kspace = tmp;
        clear tmp
    end
end

num_phase = size(downsample_kspace, 1);

%% Detect ACS lines automatically
%   Extra Nyquist-sampled k-space lines are acquired during the parallel
%   imaging scan and are used to calculate the weighting factors that
%   determine the missing k-space data.
% -------------------------------------------------------------------------

finder = logical(sum(logical(abs(downsample_kspace)), [2, 3]));
acs_finder = zeros(num_phase, 1);

% --- ACS lines are nyquist-sampled (fully-sampled k-space region), so they are adjecent.
for k = 2 : (num_phase - 1)
    if (finder(k-1) && finder(k)) || (finder(k) && finder(k+1))
        acs_finder(k) = 1;
    end
end

% --- get the indices of ACS lines
acs_lines = find(acs_finder);
num_acs_lines = numel(acs_lines);

%% Compute GRAPPA weights by solving W x A = B
%   W: weights
%   A: sources points
%   B: target points
%
%   Using ACS lines (Autocalibration signal) to calculate the weighting
%   factors that determine the missing k-space data. Note that the
%   weights are displacement-invariant.
% -------------------------------------------------------------------------

% --- size of the computing kernel
% Example: if kernel_size = [2, 3] and acceleration rate = 3, then the
% computing kernel is as below, in which the 'O' means the sources points (
% acquired data) and the 'X' represents the target points (or the
% unacquired data). Therefore, the size of this computing kernel is [4, 3].
%        O O O
%          X
%          X
%        O O O
kernel_width = kernel_size(2);
kernel_height = (kernel_size(1) - 1) * rate + 1;

assert(num_acs_lines > kernel_height, 'ACS is too short for kernel parameters')

half_width = floor(kernel_width / 2);
% half_height = floor(kernel_height / 2);

% --- parameters preparation
num_target = (rate - 1) * (kernel_size(1) - 1);                                  % number of target points in one GRAPPA computing kernel
num_source = kernel_size(1) * kernel_size(2);           % number of source points in one GRAPPA computing kernel
num_kernel_phase = num_acs_lines - kernel_height + 1;   % possible computing kernels along phase dimension
num_kernel_freq = num_freq;                             % possible computing kernels along freq dimension
num_rep = num_kernel_phase * num_kernel_freq;           % total computing kernels in ACS region

% --- define the variables of equation W x A = B
sources = zeros(num_coil * num_source, num_rep);        % matrix A
targets = zeros(num_coil * num_target, num_rep);        % matrix B

% --- Extract the ACS data and fill the matrices sources and targets
% Each computing kernel adds one column to the matrices sources and
% targets. Assume that the k-space is periodic.
count_rep = 0;
for k = acs_lines(1 : num_kernel_phase)'
    % one computing-kernel-line
    temp_source_lines = k : rate : (k + rate * (kernel_size(1) - 1));
    temp_target_lines = setdiff(k : (k + kernel_height - 1), k : rate : (k + (kernel_size(1) - 1) * rate));

    for p = 1 : num_kernel_freq
        count_rep = count_rep + 1;
        
        % Using mod() for circular indexing
        temp_source_columns = mod((p : (p + kernel_width - 1)) - 1, num_freq) + 1;
        temp_target_columns = mod(p + half_width - 1, num_freq) + 1;
        
        temp_sources = downsample_kspace(temp_source_lines, temp_source_columns, :);
        temp_targets = downsample_kspace(temp_target_lines, temp_target_columns, :);
        
        sources(:, count_rep) = reshape(temp_sources, numel(temp_sources), []);
        targets(:, count_rep) = reshape(temp_targets, numel(temp_targets), []);
    end
end

% --- compute the weights, W = B x pinv(A)
weights = targets * pinv(sources);
disp('Done£ºcompute the weights.');

%% GRAPPA reconstruction
% --- find the lines to be filled
% Corresponds to the first line to be filled of each computing kernel 
fill_finder = find(diff(finder) == -1) + 1;

recover_kspace = zeros(size(downsample_kspace));

for k = fill_finder'
    temp_source_lines = mod(((k-1) : rate : (k-1 + rate * (kernel_size(1) - 1))) - 1, num_phase) + 1;
    % lines to be filled
    temp_target_lines = setdiff((k-1) : (k + kernel_height - 2), (k-1) : rate : (k-1 + (kernel_size(1) - 1) * rate));
    for p = 1 : num_kernel_freq
        temp_source_columns = mod((p : (p + kernel_width - 1)) - 1, num_freq) + 1;
        temp_target_columns = mod(p + half_width - 1, num_freq) + 1;
        
        temp_sources = downsample_kspace(temp_source_lines, temp_source_columns, :);
        temp_targets = weights * reshape(temp_sources, numel(temp_sources), []);
        
        temp_targets = reshape(temp_targets, (rate-1)*(kernel_size(1)-1), 1, num_coil);
        recover_kspace(temp_target_lines, temp_target_columns, :) = temp_targets;
    end
end

recover_kspace(finder, :, :) = downsample_kspace(finder, :, :);

%% return
recover_kspace = recover_kspace(1:num_phase_orig, :, :);
recover_kspace = autopermute(recover_kspace, dim);
disp('Congradulation! GRAPPA reconstruction done.')

end