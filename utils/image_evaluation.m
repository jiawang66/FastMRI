function [MSE, SNR, PSNR, SSIM] = image_evaluation(x, y)

x = mat2gray(abs(x)) * 255;
y = mat2gray(abs(y)) * 255;

A = y - x;
B = x .* y;

% --- MSE, SNR, PSNR
MSE = sum(A(:).*A(:)) / numel(y);
SNR = 10 * log10(sum(x(:).*x(:)) / MSE / numel(y));
PSNR = 10 * log10(255^2 / MSE);

% --- SSIM
% ux = sum(x(:).*x(:)) / numel(x);
% uy = sum(y(:).*y(:)) / numel(y);
% sigmodx = sum(x(:).*x(:) - ux) / numel(x);
% sigmody = sum(y(:).*y(:) - uy) / numel(y);
% sigmodxy = sum(B(:).*B(:)) / (numel(B)*ux*uy) - ux*uy;
% SSIM = (2*ux*uy)*(2*sigmodxy)/(ux*ux+uy*uy)/(sigmodx*sigmodx+sigmody*sigmody);
SSIM = myssim(x,y);

end