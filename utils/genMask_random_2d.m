function mask_var = genMask_random_2d(pdf, shape, rate, hw)

randpe = round(anyrand(round(shape(1)/rate), pdf) * shape(1));
randro = round(anyrand(round(shape(2)/rate), pdf) * shape(2));
mask_var = zeros(shape);

for k = 1 : numel(randpe)
    mask_var(randpe(k), :) = 1;
end

for k = 1 : numel(randro)
    mask_var(:, randro(k)) = 1;
end

hFOV = round(shape(1)/2);
mask_var((hFOV-hw) : (hFOV+hw), :) = 1;
hFOV = round(shape(2)/2);
mask_var(:, (hFOV-hw) : (hFOV+hw)) = 1;

end