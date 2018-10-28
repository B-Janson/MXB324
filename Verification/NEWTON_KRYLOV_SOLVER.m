close all
clear
clc
format compact
%% Part 0 Initialisation
% Get the parameters to solve for
[PARAMS_SSL] = VERIF_PARAMS;
save('PARAMS_SSL')
% Get the grid & other info about grid
[DIM_SSL] = VERIF_GRID_COORD(PARAMS_SSL);
save('DIM_SSL')
% Get initial conditions
[h_old, S_old, phi_old, k_old] = VERIF_INIT_COND(DIM_SSL);

h = h_old;
h_store_SSL=h_old;
% Initial calculation of F at t = dt
F = VERIF_FVM(DIM_SSL, h, h_old, S_old, phi_old, k_old, PARAMS_SSL.dt, PARAMS_SSL);
err = norm(F, 2);
err_old = err;
PARAMS_SSL.F0=norm(F,inf);
framenum = 1;
F_old=F;
% Get the jacobian
J = JAC_FUNC(DIM_SSL, F, @VERIF_FVM, h, h_old, S_old, phi_old, k_old, PARAMS_SSL.dt, PARAMS_SSL);
M = J;

h = h_old;
S = S_old;
phi = phi_old;
k = k_old;

T_SSL = 0;
phi_avg = PHI_AVG(DIM_SSL, phi);
phi_0 = phi_avg;
phi_true = PHI_TRUE(DIM_SSL, PARAMS_SSL, 0, phi_0);

% videoName_wcont = 'WaterContent.avi';
% videoName_phead = 'PressureHead.avi';

% wcontvideo = VideoWriter(videoName_wcont);
% pheadvideo = VideoWriter(videoName_phead);

if PARAMS_SSL.realtime_plot
    % if realtime plotting, show the initial state of solution
    figure('Position', [100 160 850 500]);
    VISUALISE(DIM_SSL, h, phi, T_SSL, phi_true, phi_avg);
end

t = 0;
timesteps = 0;
steady_state = false;
m = PARAMS_SSL.gmres_max;
Pumps_SSL=zeros(DIM_SSL.n*DIM_SSL.m,1);
Evapot_SSL=zeros(DIM_SSL.n*DIM_SSL.m,1);
%% Part 1 Main Solver
tic;
% Time step iteration
while t < PARAMS_SSL.endtime
    tic
    t = t + PARAMS_SSL.dt;
    timesteps = timesteps + 1;
    iters = 0;
    F_0_norm = norm(F, Inf);
    
    % Newton step iteration
    while err > PARAMS_SSL.tol_a + PARAMS_SSL.tol_r * err_old && iters < PARAMS_SSL.max_iters

        if m > PARAMS_SSL.gmres_max
            fprintf('Recalculating Preconditioner\n');
            M=J;
            J = JAC_FUNC(DIM_SSL, F, @VERIF_FVM, h, h_old, S_old, phi_old, k_old, t, PARAMS_SSL);
            %M=J; %dt became very small
            %[M,~] = lu(J); %Quite slow
            %[M,~] = ilu(J); %Singularity warnings
            % M=diag(J); %Invariant krylov subspace
        end
        
        % Get the del h
        [dh, m] = NEWTON_GMRES(F,F_old,iters, M, PARAMS_SSL, DIM_SSL, @VERIF_FVM, h, h_old, S_old, phi_old, k_old, t);
        
        % Update estimate for current timestep's h
        h = LineSearch(DIM_SSL, @VERIF_FVM, h, dh, h_old, S_old, phi_old, k_old, t, PARAMS_SSL);
        h_store_SSL(:,end+1)=h;
        % Update F and all other variables for this time step
        F_old = F;
        [F, S, phi, k,PUMPS, EVAPOT] = VERIF_FVM(DIM_SSL, h, h_old, S_old, phi_old, k_old, t, PARAMS_SSL);
        err = norm(F, 2);
        iters = iters + 1;
        Pumps_SSL(:,end+1)=PUMPS;
        Evapot_SSL(:,end+1)=EVAPOT;
        
        
        % If haven't converged but has been too many iterations, halve time
        % step
        if iters == PARAMS_SSL.max_iters - 1 || err > 1e12
            iters = 0;
            t = t - PARAMS_SSL.dt;
            PARAMS_SSL.dt = PARAMS_SSL.dt / 2;
            t = t + PARAMS_SSL.dt;
        end
    end
    
    % Output some debug info if wanted
    if PARAMS_SSL.debug == true
        fprintf('Newton converged with GMRES for t=%d in %d iterations. Err:%d dt:%d\n', ...
            t, iters, err, PARAMS_SSL.dt);
    end
    
    T_SSL(end + 1) = t;
    phi_avg(end + 1) = PHI_AVG(DIM_SSL, phi);
    phi_true(end + 1) = PHI_TRUE(DIM_SSL, PARAMS_SSL, t, phi_0);
    
    % We have now converged, so update variables
    h_old = h;
    S_old = S;
    phi_old = phi;
    k_old = k;

    % Recalculate base error
    F = VERIF_FVM(DIM_SSL, h, h, S, phi, k, t, PARAMS_SSL);
    err = norm(F, 2);
    err_old = err;
    
    % If adaptive time stepping and converged quickly, increase time step
    if iters <= PARAMS_SSL.gmres_max
        PARAMS_SSL.dt = min(PARAMS_SSL.dt * PARAMS_SSL.adaptive_timestep, PARAMS_SSL.max_dt);
    end
    
    if PARAMS_SSL.realtime_plot == true
        VISUALISE(DIM_SSL, h, phi, T_SSL, phi_true, phi_avg);
    end
    toc
end
disp('Steady State Reached')
toc
save('T_SSL')
save('h_store_SSL')
save('Pumps_SSL')
save('Evapot_SSL')
% CREATE_VIDEO(wcontvideo, watercontent, 20);
% CREATE_VIDEO(pheadvideo, pressurehead, 20);

