close all
clear
clc
format compact
%% Part 0 Initialisation
% Get the parameters to solve for
[PARAMS] = VERIF_PARAMS;
% Get the grid & other info about grid
[DIM] = VERIF_GRID_COORD(PARAMS);
% Get initial conditions
[h_old, S_old, phi_old, k_old] = VERIF_INIT_COND(DIM);

h = h_old;

% Initial calculation of F at t = dt
F = VERIF_FVM(DIM, h, h_old, S_old, phi_old, k_old, PARAMS.dt, PARAMS);
err = norm(F, 2);
err_old = err;
framenum = 1;

% Get the jacobian
J = JAC_FUNC(DIM, F, @VERIF_FVM, h, h_old, S_old, phi_old, k_old, PARAMS.dt, PARAMS);
M = J;

h = h_old;
S = S_old;
phi = phi_old;
k = k_old;

T = 0;
phi_avg = PHI_AVG(DIM, phi);
phi_true = PHI_TRUE(DIM, PARAMS, 0);

% videoName_wcont = 'WaterContent.avi';
% videoName_phead = 'PressureHead.avi';

% wcontvideo = VideoWriter(videoName_wcont);
% pheadvideo = VideoWriter(videoName_phead);

if PARAMS.realtime_plot
    % if realtime plotting, show the initial state of solution
    figure('Position', [100 160 850 500]);
    VISUALISE(DIM, h, phi, T, phi_true, phi_avg);
end

t = 0;
timesteps = 0;
steady_state = false;
m = PARAMS.gmres_max;

%% Part 1 Main Solver
tic;
% Time step iteration
while steady_state == false && t < PARAMS.endtime
    t = t + PARAMS.dt;
    timesteps = timesteps + 1;
    iters = 0;
    
    % Newton step iteration
    while err > PARAMS.tol_a + PARAMS.tol_r * err_old && iters < PARAMS.max_iters
        if m > PARAMS.gmres_max / 2
            fprintf('Recalculating Preconditioner\n');
            J = JAC_FUNC(DIM, F, @VERIF_FVM, h, h_old, S_old, phi_old, k_old, t, PARAMS);
            M = J;
        end
        
        % Get the del h
        [dh, m] = NEWTON_GMRES(F, M, PARAMS, DIM, @VERIF_FVM, h, h_old, S_old, phi_old, k_old, t);
        
        % Update estimate for current timestep's h
        h = LineSearch(DIM, @VERIF_FVM, h, dh, h_old, S_old, phi_old, k_old, t, PARAMS);
        
        % Update F and all other variables for this time step
        [F, S, phi, k] = VERIF_FVM(DIM, h, h_old, S_old, phi_old, k_old, t, PARAMS);
        err = norm(F, 2);
        iters = iters + 1;
        
        % If haven't converged but has been too many iterations, halve time
        % step
        if iters == PARAMS.max_iters - 1 || err > 1e12
            iters = 0;
            t = t - PARAMS.dt;
            PARAMS.dt = PARAMS.dt / 2;
            t = t + PARAMS.dt;
        end
    end
    
    % Output some debug info if wanted
    if PARAMS.debug == true
        fprintf('Newton converged with GMRES for t=%d in %d iterations. Err:%d dt:%d\n', ...
            t, iters, err, PARAMS.dt);
    end
    
    T(end + 1) = t;
    phi_avg(end + 1) = PHI_AVG(DIM, phi);
    phi_true(end + 1) = PHI_TRUE(DIM, PARAMS, t);
    
    % We have now converged, so update variables
    h_old = h;
    S_old = S;
    phi_old = phi;
    k_old = k;
    
    % Recalculate base error
    F = VERIF_FVM(DIM, h, h, S, phi, k, t, PARAMS);
    err = norm(F, 2);
    err_old = err;
    
    % If adaptive time stepping and converged quickly, increase time step
    if iters <= PARAMS.gmres_max
        PARAMS.dt = min(PARAMS.dt * PARAMS.adaptive_timestep, PARAMS.max_dt);
    end
    
    if PARAMS.realtime_plot == true
        VISUALISE(DIM, h, phi, T, phi_true, phi_avg);
    end
end
disp('Steady State Reached')
toc

% CREATE_VIDEO(wcontvideo, watercontent, 20);
% CREATE_VIDEO(pheadvideo, pressurehead, 20);

