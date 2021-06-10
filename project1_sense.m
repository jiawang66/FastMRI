% Project 1: SENSE (SENSitivity Encoding)
%   Program execution begins and ends in this file.

clc, clear, close all

%% Initialization
% --- addpath
addpath('utils')

% --- load data
load data-kspace-8channels
shape = size(full_kspace);  % [ro, pe, channel]
dims = numel(shape);

%% Compute sensitivity map
% --- reconstruction of each channel by IDFT
full_img = myifftshift(myifftn(myifftshift(full_kspace, 1:(dims-1)), 2), 1:(dims-1));
% full_img = myifftshift(idftn(myifftshift(full_kspace, 1:(dims-1)), 2), 1:(dims-1));

% --- coils combination
full_img_com = sqrt(mean(abs(full_img).^2, 3));

% --- compute sensitivity map
sensitivity_map = full_img ./ (repmat(full_img_com,[1,1,shape(3)]) + 1e-32);

%% Simulate accelerated sampled data
% --- downsample
ds_rate = 2;    % acceleration rate 
dim = 1;        % acceleration dimension
ds_kspace = kdownsample(full_kspace, ds_rate, dim);

% --- reconstruction of downsampled kspace by IDFT
ds_img = myifftshift(myifftn(myifftshift(ds_kspace, 1:(dims-1)), 2), 1:(dims-1));
% ds_img = myifftshift(idftn(myifftshift(ds_kspace, 1:(dims-1)), 2), 1:(dims-1));

% --- combine
ds_img_com = sqrt(mean(abs(ds_img).^2, 3));

%% SENSE resconstruction
[ds_img_re, gfactor] = senseKernel(ds_kspace, sensitivity_map, ds_rate, dim);

% --- compute the error
err = normabs(ds_img_re) - normabs(full_img_com);

%% display
% --- display sensitivity map
figure(), set(gcf, 'outerposition', get(0,'screensize'));
for k = 1 : shape(3)
    subplot(2, shape(3)/2, k), imshow(abs(sensitivity_map(:,:,k)), [])
    title(['Smap of channel ', num2str(k)], 'Fontsize', 20)
end
saveas(gcf, 'project1_sense_sensitivity_map.png')

% ---
figure(), set(gcf, 'outerposition', get(0,'screensize'));
subplot(231), imshow(full_img_com, []), title('Reference', 'Fontsize', 20), colorbar
subplot(232), imshow(ds_img_com, []), title(['Aliased image, R = ', num2str(ds_rate)], 'FontSize', 20), colorbar
subplot(233), imshow(abs(ds_img_re), []), title(['SENSE reconstruction, R = ', num2str(ds_rate)], 'FontSize', 20), colorbar
subplot(234), imshow(abs(gfactor), []), title(['g-factor map, R = ', num2str(ds_rate)], 'FontSize', 20), colorbar
subplot(235), imshow(abs(err), []), title(['error map, R = ', num2str(ds_rate)], 'FontSize', 20), colorbar
saveas(gcf, ['project1_sense_R-', num2str(ds_rate), '.png'])
