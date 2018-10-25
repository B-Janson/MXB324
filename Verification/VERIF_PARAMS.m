function [PARAMS] = VERIF_PARAMS()
% INIT_PARAMS Initialises the parameters of the problem and returns them in
% a struct

% Grid structure
PARAMS.uniform      = true;
PARAMS.n            = 11;
PARAMS.m            = 17;

% rainfall
PARAMS.r_f          = [0.00171,0.0171,0.000171];  % rainfall constant [normal,flood,drought]
PARAMS.r_t          = 1;        % rain type 1=normal, 2=flood, 3=drought
PARAMS.r_m          = 2;        % rain model 1=constant, 2=cosine, 3=interpol

% pumping rates
PARAMS.PUMPS        = 0;        % pumps start off
PARAMS.town_rate    = 0.25;     % rate in m^2/s to extract from town bore
PARAMS.bore_rate    = 0.125;    % rate in m^2/s to extract from crop bore

% river
PARAMS.left_river   = [60; 65; 80; 50]; % z position of river bottom and river head on LHS [0; 0] for inactive
PARAMS.right_river  = [100; 100; 80; 50];   % z position of river bottom and river head on RHS [0; 0] for inactive
PARAMS.K_r          = 0.5;      % hydraulic conductivity of the river (0 turns it off)

% time
PARAMS.dt           = 10;       % timestep size
PARAMS.max_dt       = 30;       % maximum time step size
PARAMS.endtime      = 50 * 365;    % end time
PARAMS.adaptive_timestep = 1.1; % amount to 'fast forward' time if converging quickly

% solving methods
PARAMS.method = 'column';       % how to calculate the Jacobian, full or column
PARAMS.theta        = 1;        % temporal weighting, 1-Backward Euler, 0.5-Crank-Nicholson
PARAMS.sigma        = 0;        % stream weighting, 0-Upwinding, 1-Averaging, 2-Downwinding

% tolerances
PARAMS.tol_a        = 1e-5;     % absolute tolerance of the Newton step
PARAMS.tol_r        = 1e-5;     % relative tolerance of the Newton step
PARAMS.max_iters    = 18;       % maximum number of iterations for the Newton step

% heuristics
PARAMS.rho_min      = 0.65;     % minimum rho heuristic value to re-calculate jacobian
PARAMS.jacobian_update = 5;     % how often should the Jacobian be recalculated
PARAMS.steady_state_tol = 10^-16; % tolerance between water content to determine steady state

% GMRES
PARAMS.GMRES        = true;     % whether to solve using GMRES
PARAMS.gmres_tol    = 0.5;      % Maximumn bound upon the residual
PARAMS.gmres_max    = 20;       % Maximum gmres iterations

% debug and plots
PARAMS.realtime_plot = true;    % should a plot of the solution be produced in realtime? 
PARAMS.debug = true;            % should we display some debug status info

end
