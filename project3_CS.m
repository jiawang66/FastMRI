% Project: Compressive Sensing based MR Reconstruction
%   Program execution begins and ends in this file.

clc, clear, close all

%% Initialization
% --- addpath
addpath('utils')

% --- set parameters
rate = 4;               % accelaration rate (must be an integer)
dim = 1;                % accelaration dimension, 1 or 2
acs = 24;               % number of ACS (Autocalibration Signal) line
kernel_size = [2, 3];   % size of GRAPPA kernel

%% load data
% --- load full-sampled k-space data (2-D imaging with 8 channels)
load data-kspace-8channels
[num_pe, num_ro, num_ch] = size(full_kspace);   % [readout, phase-encoding, channel]
dims = numel(size(full_kspace));                % number of dimensions

% --- reconstruct the full-sampled image
full_img = myifftshift(myifftn(myifftshift(full_kspace, 1:(dims-1)), 2), 1:(dims-1));
% --- coil combination by RSS
full_img_com = sqrt(mean(abs(full_img).^2, 3));

%% generate masks
% --- generate a mask for random cartesian downsampling
load data-PEsample-pdf
mask_cs = genMask_random_2d(pe_pdf, [num_pe, num_ro], rate, 12);

% --- generate SENSE sampling mask
mask_sense = genMask_sense_2d([num_pe, num_ro], rate, dim);

% --- generate GRAPPA sampling mask
mask_grappa = genMask_grappa_2d([num_pe, num_ro], rate, dim, acs);

%% CS-NLCG (Compressive Sensing solved by Non-Linear Conjugate Gradient)
% --- generate Fourier sampling operator
FT = p2DFT(mask_cs, [num_pe, num_ro], 1, 2);

% --- subsample the kspace
sub_kspace_random = multicoil_op_2d(FT, full_img);
% --- reconstruct the subsample image
sub_img_random = multicoil_op_2d(FT', sub_kspace_random);
% --- coil combination
sub_img_random_com = sqrt(mean(abs(sub_img_random).^2, 3));

% --- generate transform operator
wave = Wavelet('Daubechies',4,4);

% --- initialize Parameters for reconstruction
params = ncg_params_init(FT, wave);

% --- CS reconstruction by Non-linear CG
% normalize
sub_kspace_in = sub_kspace_random ./ max(abs(sub_img_random(:)));
sub_img_ncg = ncg_cs(sub_kspace_in, params);
sub_img_ncg_com = sqrt(mean(abs(sub_img_ncg).^2, 3));
% --- compute error map
sub_img_ncg_error = abs(normabs(full_img_com) - normabs(sub_img_ncg_com));

%% SENSE (SENSitivity Encoding)
% --- compute sensitivity map
sensitivity_map = genSensitivityMap_2d(full_kspace);

% --- generate Fourier sampling operator
FT = p2DFT(mask_sense, [num_pe, num_ro], 1, 2);

% --- subsample the kspace
sub_kspace_sense = multicoil_op_2d(FT, full_img);
sub_img_sense = multicoil_op_2d(FT', sub_kspace_sense);
sub_img_sense_com = sqrt(mean(abs(sub_img_sense).^2, 3));

% --- SENSE reconstruction
[sub_img_sense_re, gfactor] = senseKernel(sub_kspace_sense, sensitivity_map, rate, dim);
sub_img_sense_re_com = sqrt(mean(abs(sub_img_sense_re).^2, 3));
sub_img_sense_error = abs(normabs(full_img_com) - normabs(sub_img_sense_re_com));

%% GRAPPA (GeneRalized Autocalibrating Partially Parallel Acquisitions)
% --- generate Fourier sampling operator
FT = p2DFT(mask_grappa, [num_pe, num_ro], 1, 2);

% --- subsample the kspace
sub_kspace_grappa = multicoil_op_2d(FT, full_img);
sub_img_grappa = multicoil_op_2d(FT', sub_kspace_grappa);
sub_img_grappa_com = sqrt(mean(abs(sub_img_grappa).^2, 3));

% GRAPPA reconstruction
sub_kspace_grappa_re = grappa_2d(sub_kspace_grappa, rate, dim, kernel_size);
sub_img_grappa_re = myifftshift(myifftn(myifftshift(sub_kspace_grappa_re, 1:(dims-1)), 2), 1:(dims-1));
sub_img_grappa_re_com = sqrt(mean(abs(sub_img_grappa_re).^2, 3));
sub_img_grappa_error = abs(normabs(full_img_com) - normabs(sub_img_grappa_re_com));

%% display
figure(), set(gcf, 'outerposition', [11,11,1000,1200]);
subplot(4,4,2), imshow(mask_sense, []), title('SENSE', 'FontSize', 16)
subplot(4,4,3), imshow(mask_grappa), title('GRAPPA', 'FontSize', 16)
subplot(4,4,4), imshow(mask_cs, []), title('CS', 'FontSize', 16)
subplot(4,4,5), imshow(normabs(full_img_com))
subplot(4,4,6), imshow(normabs(sub_img_sense_com))
subplot(4,4,7), imshow(normabs(sub_img_grappa_com))
subplot(4,4,8), imshow(normabs(sub_img_random_com))
subplot(4,4,10), imshow(normabs(sub_img_sense_re_com))
subplot(4,4,11), imshow(normabs(sub_img_grappa_re_com))
subplot(4,4,12), imshow(normabs(sub_img_ncg_com))
h = subplot(4,4,14); imagesc(sub_img_sense_error); axis image; h.XAxis.Visible = 'off'; h.YAxis.Visible = 'off';
h = subplot(4,4,15); imagesc(sub_img_grappa_error); axis image; h.XAxis.Visible = 'off'; h.YAxis.Visible = 'off';
h = subplot(4,4,16); imagesc(sub_img_ncg_error); axis image; h.XAxis.Visible = 'off'; h.YAxis.Visible = 'off';
saveas(gcf, ['project_R-', num2str(rate), '_compare.png'])

%% evaluation of the reconstruction
MSE = zeros(3,1);
SNR = zeros(3,1);
PSNR = zeros(3,1);
SSIM = zeros(3,1);
[MSE(1), SNR(1), PSNR(1), SSIM(1)] = image_evaluation(full_img_com, sub_img_sense_re_com);
[MSE(2), SNR(2), PSNR(2), SSIM(2)] = image_evaluation(full_img_com, sub_img_grappa_re_com);
[MSE(3), SNR(3), PSNR(3), SSIM(3)] = image_evaluation(full_img_com, sub_img_ncg_com);

METHODS = ["SENSE", "GRAPPA", "CS"];

fid = fopen(['project_R-', num2str(rate), '_evaluation.txt'],'wt');
fprintf(fid, '| methods | MSE | SNR | PSNR | SSIM |\n');
fprintf(fid, '| -- | -- | -- | -- | -- |\n');
for k = 1 : 3
    fprintf(fid, strcat("| ", METHODS(k)));
    fprintf(fid, ' | %0.3f |  %0.3f | %0.3f | %0.3f |\n', MSE(k), SNR(k), PSNR(k), SSIM(k));
end
fclose(fid);
