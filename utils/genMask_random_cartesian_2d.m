function [mask_var, mask_pdf] = genMask_random_cartesian_2d(pdf, shape, rate, dim, hw)

randpe = round(anyrand(round(shape(1)/rate), pdf) * shape(1));
mask_var = zeros(shape);
mask_pdf = repmat(pdf, [1, shape(2)]);

for k = 1 : numel(randpe)
    mask_var(randpe(k), :) = 1;
end
hFOV = round(shape(1)/2);
mask_var((hFOV-hw) : (hFOV+hw), :) = 1;

if dim == 2
    mask_var = autopermute(mask_var, dim);
    mask_pdf = autopermute(mask_pdf, dim);
end

end