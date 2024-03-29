function [PARAMS] = VERIF_PARAMS()
% INIT_PARAMS Initialises the parameters of the problem and returns them in
% a struct
%   Detailed explanation goes here

PARAMS.r_f          = [0.00171];  % rainfall constant [normal,flood,drought]
PARAMS.r_t          = 1;        % rain type 1=normal, 2=flood, 3=drought
PARAMS.r_m          = 1;        % rain model 1=constant, 2=cosine, 3=interpol
<<<<<<< HEAD
PARAMS.dt           = 20;       % timestep size
PARAMS.max_dt       = 50;       % maximum time step size
PARAMS.endtime      = 20000;    % end time
PARAMS.theta        = 1;        % temporal weighting, 1-Backward Euler, 0.5-Crank-Nicholson
PARAMS.sigma        = 0;        % stream weighting, 0-Upwinding, 1-Averaging, 2-Downwinding
PARAMS.tol_a        = 1e-5;     % absolute tolerance of the Newton step
PARAMS.tol_r        = 0;        % relative tolerance of the Newton step
PARAMS.max_iters    = 18;       % maximum number of iterations for the Newton step
PARAMS.rho_min      = 0.65;      % minimum rho heuristic value to re-calcukate jacobian
PARAMS.jacobian_update = 5;     % how often should the Jacobian be recalculated
=======
PARAMS.dt           = 5;       % timestep size
PARAMS.max_dt       = 25;       % maximum time step size
PARAMS.endtime      = 40000;    % end time
PARAMS.theta        = 0.5;      % temporal weighting, 1-Backward Euler, 0.5-Crank-Nicholson
PARAMS.sigma        = 1;        % stream weighting, 0-Upwinding, 1-Averaging, 2-Downwinding
PARAMS.tol_a        = 1e-5;     % absolute tolerance of the Newton step
PARAMS.tol_r        = 0;        % relative tolerance of the Newton step
PARAMS.max_iters    = 18;       % maximum number of iterations for the Newton step
PARAMS.rho_max      = 0.6;
PARAMS.jacobian_update = 6;     % how often should the Jacobian be recalculated
>>>>>>> master
PARAMS.realtime_plot = true;    % should a plot of the solution be produced in realtime? 
PARAMS.gmres_tol = 10^-8;         % Maximumn bound upon the residual
PARAMS.gmres_max = 16;            % Maximum gmres iterations
PARAMS.n            = 11;
PARAMS.m            = 9;
PARAMS.town_rate    = 0.25;
PARAMS.bore_rate    = 0.125;
PARAMS.K_r          = 0.1;

PARAMS.PUMPS = 0;               % PUMPS start off
PARAMS.debug = true;            % should we display some debug status info
PARAMS.adaptive_timestep = 1.1; % amount to 'fast forward' time if converging quickly
PARAMS.method = 'column';       % how to calculate the Jacobian, Full is the normal one
PARAMS.GMRES = false;
PARAMS.steady_state_tol = 10^-8; % tolerance between water content to determine steady state

end
