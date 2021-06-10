function [c, s, y] = mydwt2(x, level, wname)
%
% function [c, s, y] = mydwt2(x, level)
%
%   Similar to dwt2(), but will do one more step: put all wavelet
%   components on a 2-D image.
%
% --- INPUT
% x: 2-D image
% level: decomposition level, < log2(N)

%% Parameters check
assert(numel(size(x)) == 2, '2-D image expected');

if nargin < 2
    level = 1;
    wname = 'haar';
elseif nargin < 3
    wname = 'haar';
end

%% Compute
% --- compute wavelet component
[c,s] = wavedec2(x, level, wname);
y = zeros(size(x));

% --- put all components together
for k = level : -1 : 1
    % --- extract the high frequency component of level k
    [h, v, d] = detcoef2('all', c, s, k);
    ind = size(s,1) - k;
    
    % --- approximate component (low frequency)
    if k == level
        a = appcoef2(c, s, wname, k);
        y(1 : s(ind,1), 1 : s(ind,2)) = a;
    end
    
    % --- horizontal high frequency component 
    y(1 : s(ind,1), (1+s(ind,2)) : s(ind,2)*2) = h;
    
    % --- vertical high frequency component 
    y((1+s(ind,1)) : s(ind,1)*2, 1 : s(ind,2)) = v;
    
    % --- diagonal high frequency component 
    y((1+s(ind,1)) : s(ind,1)*2, (1+s(ind,2)) : s(ind,2)*2) = d;
end

end
