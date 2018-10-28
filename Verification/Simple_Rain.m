%% Constant Rainfall With Pumps on
%Run from the steady state for 10 years to invesigate pumping


close all
clear
clc
format compact
%% Part 0 Initialisation
% Get the parameters to solve for
tic
 load('PARAMS_SSL')
PARAMS_SIMPLE_PUMPS=PARAMS_SSL;
PARAMS_SIMPLE_PUMPS.PUMPING=1;
PARAMS_SIMPLE_PUMPS.endtime=30*365;
PARAMS_SIMPLE_PUMPS.realtime_plot=1;
PARAMS_SIMPLE_PUMPS.PUMPS(1, 3) = 0.427;
PARAMS_SIMPLE_PUMPS.PUMPS(2, 3) = 0.213;
save('PARAMS_SIMPLE_PUMPS')
% Get the grid & other info about grid
load('DIM_SSL');
DIM_SIMPLE_PUMPS=DIM_SSL;

%Recalculate 'initial conditions' using steady state
load('h_store_SSL');
h = h_store_SSL(:,end);
h_old=h;
h_store_SIMPLE_PUMPS=h_old;
S=zeros(DIM_SIMPLE_PUMPS.n*DIM_SIMPLE_PUMPS.m,1);
phi=zeros(DIM_SIMPLE_PUMPS.n*DIM_SIMPLE_PUMPS.m,1);
k=zeros(DIM_SIMPLE_PUMPS.n*DIM_SIMPLE_PUMPS.m,1);
for i = 1:DIM_SIMPLE_PUMPS.n*DIM_SIMPLE_PUMPS.m
    S(i) = SATURATION(DIM_SIMPLE_PUMPS, h, i);
    phi(i) = WATER_CONTENT(DIM_SIMPLE_PUMPS, h, S, i);
    k(i) = PERM(DIM_SIMPLE_PUMPS, h, S, i);
end
% Initial calculation of F at t = dt
F = VERIF_FVM(DIM_SIMPLE_PUMPS, h, h_old, S_old, phi_old, k_old, PARAMS_SIMPLE_PUMPS.dt, PARAMS_SIMPLE_PUMPS);
err = norm(F, 2);
err_old = err;
PARAMS_SIMPLE_PUMPS.F0=norm(F,inf);
framenum = 1;
F_old=F;
% Get the jacobian
J = JAC_FUNC(DIM_SIMPLE_PUMPS, F, @VERIF_FVM, h, h_old, S_old, phi_old, k_old, PARAMS_SIMPLE_PUMPS.dt, PARAMS_SIMPLE_PUMPS);
M = J;

h = h_old;
S = S_old;
phi = phi_old;
k = k_old;

T_SIMPLE_PUMPS = 0;
phi_avg = PHI_AVG(DIM_SIMPLE_PUMPS, phi);
phi_0 = phi_avg;
phi_true = PHI_TRUE(DIM_SIMPLE_PUMPS, PARAMS_SIMPLE_PUMPS, 0, phi_0);

% videoName_wcont = 'WaterContent.avi';
% videoName_phead = 'PressureHead.avi';

% wcontvideo = VideoWriter(videoName_wcont);
% pheadvideo = VideoWriter(videoName_phead);

if PARAMS_SIMPLE_PUMPS.realtime_plot
    % if realtime plotting, show the initial state of solution
    figure('Position', [100 160 850 500]);
    VISUALISE(DIM_SIMPLE_PUMPS, h, phi, T_SIMPLE_PUMPS, phi_true, phi_avg);
end

t = 0;
timesteps = 0;
steady_state = false;
m = PARAMS_SIMPLE_PUMPS.gmres_max;
Pumps_SIMPLE_PUMPS=zeros(DIM_SIMPLE_PUMPS.n*DIM_SIMPLE_PUMPS.m,1);
Evapot_SIMPLE_PUMPS=zeros(DIM_SIMPLE_PUMPS.n*DIM_SIMPLE_PUMPS.m,1);
%% Part 1 Main Solver
Overhead_Time_SIMPLE_PUMPS=toc
Iteration_Time_SIMPLE_PUMPS=0;
Total_Time_SIMPLE_PUMPS=Overhead_Time_SIMPLE_PUMPS;
% Time step iteration
while t < PARAMS_SIMPLE_PUMPS.endtime
    tic
    t = t + PARAMS_SIMPLE_PUMPS.dt;
    timesteps = timesteps + 1;
    iters = 0;
    F_0_norm = norm(F, Inf);
    
    % Newton step iteration
    while err > PARAMS_SIMPLE_PUMPS.tol_a + PARAMS_SIMPLE_PUMPS.tol_r * err_old && iters < PARAMS_SIMPLE_PUMPS.max_iters

        if m > PARAMS_SIMPLE_PUMPS.gmres_max
            fprintf('Recalculating Preconditioner\n');
            M=J;
            %[M,~] = lu(J); 
            %[M,~] = ilu(J); 
            %M=diag(J); % Invariant subspace
            J = JAC_FUNC(DIM_SIMPLE_PUMPS, F, @VERIF_FVM, h, h_old, S_old, phi_old, k_old, t, PARAMS_SIMPLE_PUMPS);

        end
        
        % Get the del h
        [dh, m] = NEWTON_GMRES(F,F_old,iters, M, PARAMS_SIMPLE_PUMPS, DIM_SIMPLE_PUMPS, @VERIF_FVM, h, h_old, S_old, phi_old, k_old, t);
        
        % Update estimate for current timestep's h
        h = LineSearch(DIM_SIMPLE_PUMPS, @VERIF_FVM, h, dh, h_old, S_old, phi_old, k_old, t, PARAMS_SIMPLE_PUMPS);
        h_store_SIMPLE_PUMPS(:,end+1)=h;
        % Update F and all other variables for this time step
        F_old = F;
        [F, S, phi, k,PUMPS, EVAPOT] = VERIF_FVM(DIM_SIMPLE_PUMPS, h, h_old, S_old, phi_old, k_old, t, PARAMS_SIMPLE_PUMPS);
        err = norm(F, 2);
        iters = iters + 1;
        Pumps_SIMPLE_PUMPS(:,end+1)=PUMPS;
        Evapot_SIMPLE_PUMPS(:,end+1)=EVAPOT;
        
        
        % If haven't converged but has been too many iterations, halve time
        % step
        if iters == PARAMS_SIMPLE_PUMPS.max_iters - 1 || err > 1e12
            iters = 0;
            t = t - PARAMS_SIMPLE_PUMPS.dt;
            PARAMS_SIMPLE_PUMPS.dt = PARAMS_SIMPLE_PUMPS.dt / 2;
            t = t + PARAMS_SIMPLE_PUMPS.dt;
        end
    end
    
    % Output some debug info if wanted
    if PARAMS_SIMPLE_PUMPS.debug == true
        fprintf('Newton converged with GMRES for t=%d in %d iterations. Err:%d dt:%d\n', ...
            t, iters, err, PARAMS_SIMPLE_PUMPS.dt);
    end
    
    T_SIMPLE_PUMPS(end + 1) = t;
    phi_avg(end + 1) = PHI_AVG(DIM_SIMPLE_PUMPS, phi);
    phi_true(end + 1) = PHI_TRUE(DIM_SIMPLE_PUMPS, PARAMS_SIMPLE_PUMPS, t, phi_0);
    
    % We have now converged, so update variables
    h_old = h;
    S_old = S;
    phi_old = phi;
    k_old = k;

    % Recalculate base error
    F = VERIF_FVM(DIM_SIMPLE_PUMPS, h, h, S, phi, k, t, PARAMS_SIMPLE_PUMPS);
    err = norm(F, 2);
    err_old = err;
    
    % If adaptive time stepping and converged quickly, increase time step
    if iters <= PARAMS_SIMPLE_PUMPS.gmres_max
        PARAMS_SIMPLE_PUMPS.dt = min(PARAMS_SIMPLE_PUMPS.dt * PARAMS_SIMPLE_PUMPS.adaptive_timestep, PARAMS_SIMPLE_PUMPS.max_dt);
    end
    
    if (PARAMS_SIMPLE_PUMPS.realtime_plot == true)
        VISUALISE(DIM_SIMPLE_PUMPS, h, phi, T_SIMPLE_PUMPS, phi_true, phi_avg);
    end
    Iteration_Time_SIMPLE_PUMPS(end+1)=toc;
    disp(['Iter Time ' num2str(Iteration_Time_SIMPLE_PUMPS(end))])
    Total_Time_SIMPLE_PUMPS=Total_Time_SIMPLE_PUMPS+Iteration_Time_SIMPLE_PUMPS(end);

end

disp('Steady State Reached')
disp(['Total Time ' num2str(Total_Time_SIMPLE_PUMPS(end))])
    
save('T_SIMPLE_PUMPS')
save('h_store_SIMPLE_PUMPS')
save('Pumps_SIMPLE_PUMPS')
save('Evapot_SIMPLE_PUMPS')
save('Iteration_Time_SIMPLE_PUMPS')
save('Total_Time_SIMPLE_PUMPS')

% CREATE_VIDEO(wcontvideo, watercontent, 20);
% CREATE_VIDEO(pheadvideo, pressurehead, 20);

