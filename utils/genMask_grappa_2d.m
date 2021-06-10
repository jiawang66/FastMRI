function mask = genMask_grappa_2d(shape, rate, dim, acs)

% Generate sampling mask for GRAPPA reconstruction simulation

% --- indices of acquired lines
ind = 1 : rate : shape(dim);

% --- indices of ACS lines
start = ceil(shape(dim)/2 - acs / 2);
acs_ind = ((start + 1) : (start + acs))';

% --- make a mask
mask = zeros(shape);
mask(ind, :) = 1;
mask(acs_ind, :) = 1;

if dim == 2
    mask = mask';
end

end