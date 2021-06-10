% Project 2: GRAPPA
%   Program execution begins and ends in this file.

clc, clear, close all

%% Initialization
% --- addpath
addpath('utils')

% --- load data
load data-kspace-8channels
shape = size(full_kspace);  % [ro, pe, channel]
dims = numel(shape);

% --- set the parameters
% simulate downsample
rate = 2;    % acceleration rate 
dim = 1;     % acceleration dimension

% number of ACS lines
num_acs_line = 24;

% GRAPPA kernel
kernel_size = [2, 3];

%% Simulate accelerated sampled data
% [ds_kspace, full_kspace] = grappa_kdownsample(full_kspace, rate, dim);
ds_kspace = kdownsample(full_kspace, rate, dim);

full_img = myifftshift(myifftn(myifftshift(full_kspace, 1:(dims-1)), 2), 1:(dims-1));
full_img_com = sqrt(mean(abs(full_img).^2, 3));

% --- ACS lines
shape = size(full_kspace);
start = ceil(shape(dim)/2 - num_acs_line / 2);
acs_line = (start + 1 : (start + num_acs_line))';

ds_kspace = autopermute(ds_kspace, dim);
full_kspace = autopermute(full_kspace, dim);
ds_kspace(acs_line, :, :) = full_kspace(acs_line, :, :);
ds_kspace = autopermute(ds_kspace, dim);
full_kspace = autopermute(full_kspace, dim);

% --- reconstruction of downsampled kspace by IDFT
ds_img = myifftshift(myifftn(myifftshift(ds_kspace, 1:(dims-1)), 2), 1:(dims-1));

% --- combine
ds_img_com = sqrt(mean(abs(ds_img).^2, 3));

%% GRAPPA reconstruction
recover_kspace = grappa_2d(ds_kspace, rate, dim, kernel_size);
recover_img = myifftshift(myifftn(myifftshift(recover_kspace, 1:(dims-1)), 2), 1:(dims-1));
recover_img_com = sqrt(mean(abs(recover_img).^2, 3));

%% display
figure(), set(gcf, 'outerposition', get(0,'screensize'));
for k = 1 : shape(3)
    subplot(2, shape(3)/2, k), imshow(normabs(ds_img(:,:,k)), []), colorbar
    title(['aliased image: channel-', num2str(k)], 'FontSize', 16)
end
saveas(gcf, ['project2_grappa_R-', num2str(rate), '_ACSLine-', num2str(num_acs_line), '_aliased_images.png'])

figure(), set(gcf, 'outerposition', get(0,'screensize'));
for k = 1 : shape(3)
    subplot(2, shape(3)/2, k), imshow(normabs(recover_img(:,:,k)), []), colorbar
    title(['GRAPPA reconstruction: channel-', num2str(k)], 'FontSize', 16)
end
saveas(gcf, ['project2_grappa_R-', num2str(rate), '_ACSLine-', num2str(num_acs_line), '_grappa_reconstruction.png'])

error = normabs(full_img_com) - normabs(recover_img_com);

figure(), set(gcf, 'outerposition', [1, 1, 1350, 1000]);
subplot(221), imshow(normabs(full_img_com), []), title('Full image', 'FontSize', 16), colorbar
subplot(222), imshow(normabs(ds_img_com), []), title('Downsampled image', 'FontSize', 16), colorbar
subplot(223), imshow(normabs(recover_img_com), []), title(['GRAPPA reconstruction: R = ', num2str(rate), ', ACSLine = ', num2str(num_acs_line)], 'FontSize', 16), colorbar
subplot(224), imshow(error, [0, 0.3]), title('Error map', 'FontSize', 16), colorbar
saveas(gcf, ['project2_grappa_R-', num2str(rate), '_ACSLine-', num2str(num_acs_line), '_comparison.png'])
