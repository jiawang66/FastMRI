function y = multicoil_op_2d(op, x, pdf)

if nargin < 3
    pdf = ones(size(x,1), size(x,2));
end

y = zeros(size(x));
for k = 1 : size(x,3)
    y(:,:,k) = op * (x(:,:,k) ./ pdf);
end

end