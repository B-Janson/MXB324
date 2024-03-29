function [PARAMS] = INIT_PARAMS()
% INIT_PARAMS Initialises the parameters of the problem and returns them in
% a struct
%   Detailed explanation goes here

PARAMS.r_f          = [0.00171,0.0171,0.000171];  % rainfall constant [normal,flood,drought]
PARAMS.r_t          = 1;        % rain type 1=normal, 2=flood, 3=drought
PARAMS.r_m          = 1;        % rain model 1=constant, 2=cosine, 3=interpol
PARAMS.K_r          = 0.1;      % river hydraulic conductivity
PARAMS.dt           = 10;       % timestep size
PARAMS.max_dt       = 30;       % maximum time step size
PARAMS.endtime      = 20000;    % end time
PARAMS.theta        = 0.5;      % temporal weighting, 1-Backward Euler, 0.5-Crank-Nicholson
PARAMS.sigma        = 1;        % stream weighting, 0-Upwinding, 1-Averaging, 2-Downwinding
PARAMS.tol_a        = 1e-10;    % absolute tolerance of the Newton step
PARAMS.tol_r        = 0;        % relative tolerance of the Newton step
PARAMS.max_iters    = 18;       % maximum number of iterations for the Newton step
PARAMS.rho_max      = 0.6;
PARAMS.jacobian_update = 5;     % how often should the Jacobian be recalculated
PARAMS.plot_times = 0:5:PARAMS.endtime; % vector of times to save and plot the numerical solution
PARAMS.realtime_plot = true;    % should a plot of the solution be produced in realtime? 
PARAMS.gmres_tol = 1e-8;         % Maximumn bound upon the residual
PARAMS.gmres_max = 16;            % Maximum gmres iterations

PARAMS.PUMPS=0;                 % PUMPS start off    %No, way to slow
PARAMS.town_rate    = 0.25;     % rate for town pump (m^2/day)
PARAMS.bore_rate    = 0.125;    % rate for crop pump (m^2/day)
PARAMS.debug = true;            % should we display some debug status info
PARAMS.adaptive_timestep = 1.2; % should we change the time step
PARAMS.method = 'column';       % how to calculate the Jacobian, full is the normal one, column is tridiagonal column
PARAMS.steady_state_tol = 1e-8; % tolerance between water content to determine steady state


end
