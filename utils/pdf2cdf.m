function cdf = pdf2cdf(pdf)
%
% function cdf = pdf2cdf(pdf)
%

pdf = pdf / sum(pdf, 'all');
cdf = zeros(size(pdf));
cdf(1) = pdf(1);

for k = 2 : numel(pdf)
    cdf(k) = cdf(k-1) + pdf(k);
end

end