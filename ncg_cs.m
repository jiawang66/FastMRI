function res = ncg_cs(kspace, params)

% Use Nonlinear Conjugate Gradient Method to solve the CS model of MRI
% reconstruction

fprintf('\n================ Start CS reconstruction ================\n');

res = zeros(size(kspace));
for k = 1 : size(kspace,3)
    % --- scale data
    im = params.FT' * kspace(:,:,k);
    params.data = kspace(:,:,k);
    params.lambda = params.lambda * max(abs(im(:)));
    im_wave = params.wave * im;
    
    figure(101), set(gcf, 'outerposition', [800,60,500,500])
    imshow(abs(im),[]), title(['channel-', num2str(k), ' initial']), drawnow
    
    % --- compute, N-step restart
    for n = 1 : params.Nsteps
        fprintf('---------- channel-%d , loop-%d ----------\n', k, n);
        im_wave = ncg_cs_kernel(im_wave, params);
        im = params.wave' * im_wave;
        figure(101), set(gcf, 'outerposition', [800,60,500,500]),
        imshow(abs(im),[]), title(['channel-', num2str(k), ', loop-', num2str(n)]), drawnow
    end
    res(:,:,k) = params.wave' * im_wave;
end

fprintf('================ CS reconstruction Done ================\n');

end

function x0 = ncg_cs_kernel(x, params)
%
% function m = ncg_cs(x0, params)
%
%   Nonlinear Conjugate Gradient Method
%
% f(x) = ||F * W' * x - y||^2 + lambda * ||x||_1 (x in wavelet domain)
% df/dx = 2 * W * F' * (F * W' * x - y) + lambda * (x' * x + mu)^(-1/2) * x

x0 = x;

% --- line search parameters
alpha = params.lineSearchAlpha;
beta = params.lineSearchBeta;
t0 = params.lineSearchT0;

% --- compute first search direction
g0 = gradientObjective(x0, params);
dx = -g0;

% --- iteration
k = 0;
while (k < params.iters) && (norm(g0(:)) > params.gradToll)
    %% backtracking linesearch   
    c = 0;
    t = t0;
    x1 = x0 + t * dx;
    f0 = objective(x0, params);
    f1 = objective(x1, params);
    while (c < params.lineSearchIters) && (f1 > f0 - alpha * t * abs(g0(:)'*dx(:)))
        c = c + 1;
        t = t * beta;
        x1 = x0 + t * dx;
        f1 = objective(x1, params);
    end
    
    D = f0 - f1;
    fprintf('Iteration-%d,\t descend = %0.4f,\t abs(grad) = %0.4f \n', k+1, D, norm(g0(:)))
    
    %% conjugate gradient calculation
    g1 = gradientObjective(x1, params);
    bk = g1(:)' * g1(:) / (g0(:)' * g0(:) + eps);
    g0 = g1;
    x0 = x1;
    dx = -g1 + bk * dx;
    k = k + 1;
end

end

%% compute the objective function
function res = objective(x, params)

% compute the objective function
% f(x) = ||F * W' * x - y||^2 + lambda * ||x||_1

p = params.pnorm;
obj = params.FT * (params.wave' * x) - params.data;
obj = obj(:)' * obj(:);

if params.lambda
    wave = (x .* conj(x) + params.smooth).^(p/2);
else
    wave = 0;
end

wave = sum(wave(:)) * params.lambda;

res = obj + wave;

end


%% compute the gradient of the objective function
function grad = gradientObjective(x, params)

% compute the gradient of the objective function
% grad = "gradient of the data consistency" + lambda * "gradient of the L1 transform operator"

gradObj = gradientObj(x, params);
gradWave = 0;

if params.lambda
    gradWave = gradientWave(x, params);
end

grad = gradObj + params.lambda * gradWave;

end

%% compute the gradient of the data consistency
function gradObj = gradientObj(x, params)

% compute the gradient of the data consistency
% gradObj = 2 * W * F' * (F * W' * x - y)

gradObj = params.FT * (params.wave' * x) - params.data;
gradObj = 2 * (params.wave * (params.FT' * gradObj));

end

%% compute the gradient of the L1 transform operator
function gradWave = gradientWave(x, params)

% compute the gradient of the L1 transform operator
% gradWave = (x' * x + mu)^(-1/2) * x

p = params.pnorm;

gradWave = p * x .* (x .* conj(x) + params.smooth).^(p/2-1);

end
