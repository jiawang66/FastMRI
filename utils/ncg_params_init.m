function params = ncg_params_init(FT, wave)

% Initial the parameters for the iteration of Non-linear Conjugate Gradient

params.FT = FT;
params.wave = wave;
params.data = 0;
params.lambda = 0.01;
params.iters = 8;
params.Nsteps = 5;
params.gradToll = 1e-4;
params.smooth = 1e-15;
params.pnorm = 1;

% --- linesearch parameters
params.lineSearchAlpha = 0.05;
params.lineSearchBeta = 0.6;
params.lineSearchT0 = 1;
params.lineSearchIters = 150;

end