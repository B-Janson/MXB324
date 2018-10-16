function [PARAMS] = INIT_PARAMS()
%INIT_PARAMS Initialises the parameters of the problem and returns them in
%a struct
%   Detailed explanation goes here

PARAMS.r_f          = [0.00171,0.0171,0.000171];  % rainfall constant [normal,flood,drought]
PARAMS.r_t          = 1;        % 1=constant, 2=cosine, 3=interpol
PARAMS.dt           = 10;       % timestep size
PARAMS.endtime      = 4000000;  % end time
PARAMS.theta        = 0.5;      % temporal weighting, 1-Backward Euler, 0.5-Crank-Nicholson
PARAMS.sigma        = 1;        % stream weighting, 0-Upwinding, 1-Averaging, 2-Downwinding
PARAMS.tol_a        = 1e-10;    % absolute tolerance of the Newton step
PARAMS.tol_r        = 0;        % relative tolerance of the Newton step
PARAMS.max_iters    = 20;       % maximum number of iterations for the Newton step
PARAMS.jacobian_update = 6;     % how often should the Jacobian be recalculated
PARAMS.plot_times = 0:5:PARAMS.endtime; % vector of times to save and plot the numerical solution
PARAMS.realtime_plot = true;    % should a plot of the solution be produced in realtime? 
PARAMS.gmres_tol=10^-8;         % Maximumn bound upon the residual
PARAMS.gmres_max=16;            % Maximum gmres iterations

PARAMS.PUMPS=0;%PUMPS start off                                %No, way to slow
PARAMS.debug = true;            % should we display some debug status info
PARAMS.method = 'Column'; % how to calculate the Jacobian, Full is the normal one
PARAMS.breaktol=10^-12;



end

