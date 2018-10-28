%Nonlinear mesh
close all
clear
clc
format compact
%% Part 0 Initialisation
% Get the parameters to solve for
[PARAMS_NL] = VERIF_PARAMS;
PARAMS_NL.uniform=0;
save('PARAMS_NL')
% Get the grid & other info about grid
[DIM_NL] = VERIF_GRID_COORD(PARAMS_NL);
save('DIM_NL')
% Get initial conditions
[h_old, S_old, phi_old, k_old] = VERIF_INIT_COND(DIM_NL);

h = h_old;
h_store_NL=h_old;
% Initial calculation of F at t = dt
F = VERIF_FVM(DIM_NL, h, h_old, S_old, phi_old, k_old, PARAMS_NL.dt, PARAMS_NL);
err = norm(F, 2);
err_old = err;
PARAMS_NL.F0=norm(F,inf);
framenum = 1;
F_old=F;
% Get the jacobian
J = JAC_FUNC(DIM_NL, F, @VERIF_FVM, h, h_old, S_old, phi_old, k_old, PARAMS_NL.dt, PARAMS_NL);
M = J;

h = h_old;
S = S_old;
phi = phi_old;
k = k_old;

T_NL = 0;
phi_avg = PHI_AVG(DIM_NL, phi);
phi_0 = phi_avg;
phi_true = PHI_TRUE(DIM_NL, PARAMS_NL, 0, phi_0);

% videoName_wcont = 'WaterContent.avi';
% videoName_phead = 'PressureHead.avi';

% wcontvideo = VideoWriter(videoName_wcont);
% pheadvideo = VideoWriter(videoName_phead);

if PARAMS_NL.realtime_plot
    % if realtime plotting, show the initial state of solution
    figure('Position', [100 160 850 500]);
    VISUALISE(DIM_NL, h, phi, T_NL, phi_true, phi_avg);
end

t = 0;
timesteps = 0;
steady_state = false;
m = PARAMS_NL.gmres_max;
Pumps_NL=zeros(DIM_NL.n*DIM_NL.m,1);
Evapot_NL=zeros(DIM_NL.n*DIM_NL.m,1);
%% Part 1 Main Solver
Overhead_Time_NL=toc
Iteration_Time_NL=0;
Total_Time_NL=Overhead_Time_NL;
% Time step iteration
while t < PARAMS_NL.endtime
    tic
    t = t + PARAMS_NL.dt;
    timesteps = timesteps + 1;
    iters = 0;
    F_0_norm = norm(F, Inf);
    
    % Newton step iteration
    while err > PARAMS_NL.tol_a + PARAMS_NL.tol_r * err_old && iters < PARAMS_NL.max_iters

        if m > PARAMS_NL.gmres_max
            fprintf('Recalculating Preconditioner\n');
            M=J;
            %[M,~] = lu(J); 
            %[M,~] = ilu(J); 
            %M=diag(J); % Invariant subspace
            J = JAC_FUNC(DIM_NL, F, @VERIF_FVM, h, h_old, S_old, phi_old, k_old, t, PARAMS_NL);

        end
        
        % Get the del h
        [dh, m] = NEWTON_GMRES(F,F_old,iters, M, PARAMS_NL, DIM_NL, @VERIF_FVM, h, h_old, S_old, phi_old, k_old, t);
        
        % Update estimate for current timestep's h
        h = LineSearch(DIM_NL, @VERIF_FVM, h, dh, h_old, S_old, phi_old, k_old, t, PARAMS_NL);
        h_store_NL(:,end+1)=h;
        % Update F and all other variables for this time step
        F_old = F;
        [F, S, phi, k,PUMPS, EVAPOT] = VERIF_FVM(DIM_NL, h, h_old, S_old, phi_old, k_old, t, PARAMS_NL);
        err = norm(F, 2);
        iters = iters + 1;
        Pumps_NL(:,end+1)=PUMPS;
        Evapot_NL(:,end+1)=EVAPOT;
        
        
        % If haven't converged but has been too many iterations, halve time
        % step
        if iters == PARAMS_NL.max_iters - 1 || err > 1e12
            iters = 0;
            t = t - PARAMS_NL.dt;
            PARAMS_NL.dt = PARAMS_NL.dt / 2;
            t = t + PARAMS_NL.dt;
        end
    end
    
    % Output some debug info if wanted
    if PARAMS_NL.debug == true
        fprintf('Newton converged with GMRES for t=%d in %d iterations. Err:%d dt:%d\n', ...
            t, iters, err, PARAMS_NL.dt);
    end
    
    T_NL(end + 1) = t;
    phi_avg(end + 1) = PHI_AVG(DIM_NL, phi);
    phi_true(end + 1) = PHI_TRUE(DIM_NL, PARAMS_NL, t, phi_0);
    
    % We have now converged, so update variables
    h_old = h;
    S_old = S;
    phi_old = phi;
    k_old = k;

    % Recalculate base error
    F = VERIF_FVM(DIM_NL, h, h, S, phi, k, t, PARAMS_NL);
    err = norm(F, 2);
    err_old = err;
    
    % If adaptive time stepping and converged quickly, increase time step
    if iters <= PARAMS_NL.gmres_max
        PARAMS_NL.dt = min(PARAMS_NL.dt * PARAMS_NL.adaptive_timestep, PARAMS_NL.max_dt);
    end
    
    if (PARAMS_NL.realtime_plot == true)
        VISUALISE(DIM_NL, h, phi, T_NL, phi_true, phi_avg);
    end
    Iteration_Time_NL(end+1)=toc;
    disp(['Iter Time ' num2str(Iteration_Time_NL(end))])
    Total_Time_NL=Total_Time_NL+Iteration_Time_NL(end);

end

disp('Steady State Reached')
disp(['Total Time ' num2str(Total_Time_NL(end))])
    
save('T_NL')
save('h_store_NL')
save('Pumps_NL')
save('Evapot_NL')
save('Iteration_Time_NL')
save('Total_Time_NL')

% CREATE_VIDEO(wcontvideo, watercontent, 20);
% CREATE_VIDEO(pheadvideo, pressurehead, 20);

