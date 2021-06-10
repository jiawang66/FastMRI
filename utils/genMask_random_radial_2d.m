function mask_var = genMask_random_radial_2d(pdf, shape, rate, hw)

mask_var= zeros(shape);
N = sqrt(shape(:)'*shape(:));
num = 402;  % number of spokes for full-sampling
num = ceil(num / rate);
theta = rand(num,1) * pi;
for k = 1 : num
    index = round(anyrand(round(N), pdf) * N - N/2);
    index = union(index, -hw:hw);
    x = round(index * sin(theta(k)) + shape(1)/2);
    y = round(index * cos(theta(k)) + shape(2)/2);
    
    for s = 1 : numel(x)
        if 1 <= x(s) && x(s) <= shape(1) && 1 <= y(s) && y(s) <= shape(2)
            mask_var(x(s),y(s)) = 1;
        end
    end
end

end