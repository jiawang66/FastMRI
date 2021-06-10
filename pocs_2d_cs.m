function re = pocs_2d_cs(kspace, mask_var, mask_pdf, lambda, iter, epsilon)
%
% function re = pocs_2d_cs(x, lambda, iter)
%
% kspace: subsampled kspace with zero padding at non-acquisition points 
%

%%

num_ch = size(kspace, 3);
re = zeros(size(kspace));

wave = Wavelet('Daubechies',4,4);

for ch = 1 : num_ch
    DATA = kspace(:, :, ch) .* mask_var;
    im_tmp = ifft2c(DATA ./ mask_pdf);  % initial value
    for k = 1 : iter
        im_cs = wave * im_tmp;
        im_cs = SoftThresh(im_cs, lambda);
        
        tmp = im_tmp;
        im_tmp = wave' * im_cs;
        im_tmp = ifft2c(fft2c(im_tmp) .* (1 - mask_var) + DATA);
        error = matnorm(normabs(tmp) - normabs(im_tmp), 3);
        fprintf('lambda-%0.2f, channel-%d, iteration-%d, abs(x_{k+1} - x_{k}) = %0.4f\n', lambda, ch, k, error);
        
        if error < epsilon
            break
        end
    end
    re(:,:,ch) = im_tmp;
end

re(isnan(re)) = 0;

end