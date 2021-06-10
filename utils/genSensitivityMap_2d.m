function smap = genSensitivityMap_2d(kspace, sLine)

% simulate the acquisition of the sensitivity map

shape = size(kspace);
dims = numel(shape);

if nargin < 2
    sLine = ceil(shape(1)/5);
end

% --- indices of acquisition lines for computing sensitivity map
start = ceil(shape(1)/2 - sLine/2);
indices = ((start + 1) : (start + sLine))';

smap_kspace = zeros(shape);
smap_kspace(indices,:,:) = kspace(indices,:,:);

smap_img = myifftshift(myifftn(myifftshift(smap_kspace, 1:(dims-1)), 2), 1:(dims-1));
smap_img_com = sqrt(mean(abs(smap_img).^2, 3));
smap = smap_img ./ (repmat(smap_img_com,[1,1,shape(3)]) + eps);

end