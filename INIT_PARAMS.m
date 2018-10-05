function [PARAMS] = INIT_PARAMS()
%INIT_PARAMS Initialises the parameters of the problem and returns them in
%a struct
%   Detailed explanation goes here

PARAMS.r_f = 0.00171;         % rainfall constant
PARAMS.dt = 50;        % timestep size
PARAMS.endtime = 4000000;   % end time
PARAMS.theta = 0.5;     % temporal weighting, 1-Backward Euler, 0.5-Crank-Nicholson
PARAMS.sigma = 1;       % stream weighting, 0-Upwinding, 1-Averaging, 2-Downwinding
PARAMS.tol_a = 1e-10;   % absolute tolerance of the Newton step
PARAMS.tol_r = 0;       % relative tolerance of the Newton step
PARAMS.T = PARAMS.dt:PARAMS.dt:PARAMS.endtime; % each time step to solve for
PARAMS.max_iters = 20; % maximum number of iterations for the Newton step
PARAMS.jacobian_update = 1; % how often should the Jacobian be recalculated
PARAMS.plot_times = 0:5:PARAMS.endtime; % vector of times to save and plot the numerical solution
PARAMS.realtime_plot = true; % should a plot of the solution be produced in realtime?
PARAMS.debug = true; % should we display some debug status info
PARAMS.jacobian_method = 'finite difference'; % how to calculate the Jacobian
% 'finite differences' - Uses full finite difference jacobian requiring N function evaluations. 
% 'tridiagonal' - Uses a banded tridiagonal Jacobian with 3 function evaluations.

end

