function [re, gfactor] = senseKernel(dskspace, sm, r, dim, phi)
% re = senseKernel(dskspace, sm, r, dim, phi)
%   SENSE (SENSitivity Encoding) reconstruction algorithm of MRI
%
% - Input
% dskspace: aliasing kspace, 3D, 4D or 5D matrix, the last dimension should
%           be channels (coils)
% sm: sensitivity map of each coil
% r: acceleration rate, 1 <= r <= number of channels
% dim: acceleration dimension
% phi: coil noise correction matrix
%
% - Return
% re: reconstructed image
% gfactor: g-factor map

%% Parameters preparation
if size(dskspace) ~= size(sm)
    error('Shape of matrix kdata and sensitivity maps expected to be equal') 
end

shape = size(dskspace);

if dim > numel(shape) - 1
    error('Wrong acceleration dimension, which cannot be the last dimension of first parameter');
end

if nargin < 5
    phi = eye(shape(end));
end

%% Fourier reconstruction of downsampled kspace
dims = numel(size(dskspace));
img = myifftshift(myifftn(myifftshift(dskspace, 1:(dims-1)), numel(shape) - 1), 1:(dims-1));
% img = myifftshift(idftn(myifftshift(dskspace, 1:(dims-1)), numel(shape) - 1), 1:(dims-1));

% --- reshape the img and sm,
%   [acceleration dimenison, [], channels]
img = autopermute(img, dim);
sm = autopermute(sm, dim);
shape = size(img);

if numel(shape) > 3
    img = reshape(img, shape(1), numel(img)/shape(1)/shape(end), []);
    sm = reshape(sm, shape(1), numel(sm)/shape(1)/shape(end), []);
end

%% reconstruction
shape_re = size(img);
rows = shape_re(1) / r;         % period of aliasing image in accelerated dimension

re = zeros(shape_re(1:end-1));  % unfolded image
gfactor = zeros(shape_re(1:end-1));
index = zeros(r, 1);

flag = zeros(shape_re(1), 1);

tic
for k = 1 : rows
    % --- compute the index of unfolded voxels
    for i = 1 : r
        index(i) = round(k + (i - 1) * rows);
        flag(index(i)) = 1;
    end
    
    C = permute(sm(index, :, :), [3, 1, 2]);
    I = permute(squeeze(img(k, :, :)), [2, 1]);
    
    for m = 1 : shape_re(2)
        tempC_pinv = pinv(C(:,:,m)' * pinv(phi) * C(:,:,m));
        re(index, m) = tempC_pinv * C(:,:,m)' * pinv(phi) * I(:,m);     
        gfactor(index, m) = sqrt(diag(tempC_pinv) .* diag(C(:,:,m)' * pinv(phi) * C(:,:,m)));
    end
end
toc

% --- Interpolation
%   When acceration rate R is odd, some line may be empty after SENSE reconstruction
if sum(flag) ~= numel(flag)
    [re, gfactor] = interpolation(re, gfactor, flag, rows);
end

% --- recover the shape
re = autopermute(reshape(re, shape(1:end-1)), dim);
gfactor = autopermute(reshape(gfactor, shape(1:end-1)), dim);

end


%% Interpolate one line
function [re, gfactor] = interpolation(re, gfactor, flag, rows)

for k = ceil(rows) : numel(flag)
    if flag(k) == 0
        n = k - 1;  % index of lower nearest line be unfolded
        for m = k + 1 : numel(flag) % index of upper nearest line be unfolded
            if flag(m) == 1
                break;
            end
        end
        if isempty(m)
            m = n;
        end

        % --- Interpolation
        if m == n
            re(k, :) = re(n, :);
            gfactor(k, :) = gfactor(n, :);
            flag(k) = 1;
        else
            ind = 1 : (m - n + 1);

            tmp = interpolationKernel([1, numel(ind)]', re([n,m],:), ind, 'pchip');
            re(k:(m-1),:) = tmp(2:(end-1), :);

            tmp = interpolationKernel([1, numel(ind)]', gfactor([n,m],:), ind, 'pchip');
            gfactor(k:(m-1),:) = tmp(2:(end-1), :);

            flag(k:(m-1),:) = 1;
        end
    end
end
end


%% Interpolation complex
function y = interpolationKernel(x, v, xq, method)

y = interp1(x, v, xq, method);
y_abs = interp1(x, abs(v), xq, method);
y = y_abs .* y ./ abs(y + 1e-31);

end