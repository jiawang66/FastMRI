function y = anyrand(n, pdf)
%
% function y = anyrand(n, pdf)
%

x = rand(n, 1); % uniformly distributed pseudo-random numbers in [0,1]
y = zeros(n, 1);
cdf = pdf2cdf(pdf);

for k = 1 : n
    if x(k) <= cdf(1)
        y(k) = 1;
        continue
    elseif x(k) >= cdf(end)
        y(k) = numel(cdf);
        continue
    end
    
    for w = 1 : numel(cdf)-1
        if x(k) >= cdf(w) && x(k) < cdf(w+1)
            y(k) = w;
            break
        end
    end
end

y = y / numel(cdf);

end