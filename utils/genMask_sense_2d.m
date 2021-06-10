function mask = genMask_sense_2d(shape, rate, dim)

% Generate sampling mask for SENSE reconstruction simulation

% --- indices of acquired lines
ind = 1 : rate : shape(dim);

% --- make a mask
mask = zeros(shape);
mask(ind, :) = 1;

if dim == 2
    mask = mask';
end

end