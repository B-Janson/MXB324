close all
clear
clc
format compact
%% Part 0 Initialisation
% Get the parameters to solve for
[PARAMS] = VERIF_PARAMS;
% Get the boundary conditions
[BC] = BOUNDARY_CONDITIONS(PARAMS);
% Get the grid & other info about grid
[DIM] = VERIF_GRID_COORD(PARAMS, BC);
% Get initial conditions
[h_old, S_old, phi_old, k_old] = VERIF_INIT_COND(DIM);

h = h_old;

% Initial calculation of F at t = dt
F = VERIF_FVM(DIM, h, h_old, S_old, phi_old, k_old, PARAMS.dt, PARAMS);
err = norm(F, 2);
err_old = err;

iters = 0;
fevals = 0;
f_eval_total = 1;
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

sat_col = SAT_COLOUR;

if PARAMS.realtime_plot
    % if realtime plotting, show the initial state of solution
    figure('Position', [100 160 850 500]);
    VISUALISE(DIM, h, phi, T, phi_true, phi_avg);
end

t = 0;
timesteps = 0;
steady_state = false;

dh_guess = zeros(DIM.n * DIM.m, 1);

%% Part 1 Main Solver
tic;
while steady_state == false && t < PARAMS.endtime
    t = t + PARAMS.dt;
    timesteps = timesteps + 1;
    
    while err > PARAMS.tol_a + PARAMS.tol_r * err_old && iters < PARAMS.max_iters
        if mod(iters, PARAMS.jacobian_update) == 0 && err > PARAMS.tol_r * 10 || rho > PARAMS.rho_min
            J_old=J;
            J = JAC_FUNC(DIM, F, @VERIF_FVM, h, h_old, S_old, phi_old, k_old, t, PARAMS);
            M = ilu(J);
            fevals = fevals + DIM.n * DIM.m;
        end
        
        % Get the del h  
        if PARAMS.GMRES
            dh = NEWTON_GMRES(J, -F, dh_guess, M, PARAMS.tol_a, 20, false);
        else
            dh = J\(-F);
        end
        
        % Update estimate for current timestep's h
        h = LineSearch(DIM, @VERIF_FVM, h, dh, h_old, S_old, phi_old, k_old, t, PARAMS);

        % Update F and all other variables for this time step
        [F, S, phi, k] = VERIF_FVM(DIM, h, h_old, S_old, phi_old, k_old, t, PARAMS);
        rho = norm(F, 2) / err;
        err = norm(F, 2);
        iters = iters + 1;
        fevals = fevals + 1;
        % Output some debug info if wanted
        if PARAMS.debug == true
            fprintf('t:%d iters:%d err:%d fevals:%d timesteps:%d steady_state:%d dt:%d rho:%d norm(dh-dh2):%d\n', ...
                t, iters, err, fevals, timesteps, min(h) >= 0, PARAMS.dt, rho, norm(dh-dh,2));
        end
        
        % If haven't converged but has been too many iterations, halve time
        % step
        if iters == PARAMS.max_iters - 1 || err > 1e12
            iters = 0;
            t = t - PARAMS.dt;
            PARAMS.dt = PARAMS.dt / 2;
            t = t + PARAMS.dt;
        end
        
    end
    
    if norm(phi - phi_old, 2) < PARAMS.steady_state_tol
        fprintf('steady state\n');
        PARAMS.PUMPS = 1;
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
    
    f_eval_total = f_eval_total + fevals + 1;
    fevals = 0;
    
    % If adaptive time stepping and converged quickly, increase time step
    if iters <= PARAMS.jacobian_update
        PARAMS.dt = min(PARAMS.dt * PARAMS.adaptive_timestep, PARAMS.max_dt);
    end
    
    % reset iters
    iters = 0;

    if PARAMS.realtime_plot == true
        VISUALISE(DIM, h, phi, T, phi_true, phi_avg);
    end
    
end
disp('Steady State Reached')
toc

save('steady_state_1')

% CREATE_VIDEO(wcontvideo, watercontent, 20);
% CREATE_VIDEO(pheadvideo, pressurehead, 20);

PARAMS.PUMPS=1;

%%
clear
close all
load('steady_state_1.mat')

