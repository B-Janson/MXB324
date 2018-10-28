%% Constant Rainfall,lower conductivity With Pumps on
%Run from the steady state for 10 years to invesigate pumping


close all
clear
clc
format compact
%% Part 0 Initialisation
% Get the parameters to solve for
tic
 load('PARAMS_SSL')
PARAMS_LK=PARAMS_SSL;
PARAMS_LK.PUMPING=1;
PARAMS_LK.endtime=25*365;
PARAMS_LK.realtime_plot=1;
PARAMS_LK.K_r=0.6;
PARAMS_LK.PUMPS=[ 450, 10, 0.213; 100, 50, 0.427];
save('PARAMS_LK')
% Get the grid & other info about grid
load('DIM_SSL');
DIM_LK=DIM_SSL;

%Recalculate 'initial conditions' using steady state
load('h_store_SSL');
h = h_store_SSL(:,end);
h_old=h;
h_store_LK=h_old;
S=zeros(DIM_LK.n*DIM_LK.m,1);
phi=zeros(DIM_LK.n*DIM_LK.m,1);
k=zeros(DIM_LK.n*DIM_LK.m,1);
for i = 1:DIM_LK.n*DIM_LK.m
    S(i) = SATURATION(DIM_LK, h, i);
    phi(i) = WATER_CONTENT(DIM_LK, h, S, i);
    k(i) = PERM(DIM_LK, h, S, i);
end

% Initial calculation of F at t = dt
F = VERIF_FVM(DIM_LK, h, h_old, S_old, phi_old, k_old, PARAMS_LK.dt, PARAMS_LK);
err = norm(F, 2);
err_old = err;
PARAMS_LK.F0=norm(F,inf);
framenum = 1;
F_old=F;
% Get the jacobian
J = JAC_FUNC(DIM_LK, F, @VERIF_FVM, h, h_old, S_old, phi_old, k_old, PARAMS_LK.dt, PARAMS_LK);
M = J;

h = h_old;
S = S_old;
phi = phi_old;
k = k_old;

T_LK = 0;
phi_avg = PHI_AVG(DIM_LK, phi);
phi_0 = phi_avg;
phi_true = phi_0;

% videoName_wcont = 'WaterContent.avi';
% videoName_phead = 'PressureHead.avi';

% wcontvideo = VideoWriter(videoName_wcont);
% pheadvideo = VideoWriter(videoName_phead);

if PARAMS_LK.realtime_plot
    % if realtime plotting, show the initial state of solution
    figure('Position', [100 160 850 500]);
    VISUALISE(DIM_LK, h, phi, T_LK, phi_true, phi_avg);
end

t = 0;
timesteps = 0;
steady_state = false;
m = PARAMS_LK.gmres_max;
Pumps_LK=zeros(DIM_LK.n*DIM_LK.m,1);
Evapot_LK=zeros(DIM_LK.n*DIM_LK.m,1);
%% Part 1 Main Solver
Overhead_Time_LK=toc
Iteration_Time_LK=0;
Total_Time_LK=Overhead_Time_LK;
% Time step iteration
while t < PARAMS_LK.endtime
    tic
    t = t + PARAMS_LK.dt;
    timesteps = timesteps + 1;
    iters = 0;
    F_0_norm = norm(F, Inf);
    
    % Newton step iteration
    while err > PARAMS_LK.tol_a + PARAMS_LK.tol_r * err_old && iters < PARAMS_LK.max_iters

        if m > PARAMS_LK.gmres_max
            fprintf('Recalculating Preconditioner\n');
            M=J;
            %[M,~] = lu(J); 
            %[M,~] = ilu(J); 
            %M=diag(J); % Invariant subspace
            J = JAC_FUNC(DIM_LK, F, @VERIF_FVM, h, h_old, S_old, phi_old, k_old, t, PARAMS_LK);

        end
        
        % Get the del h
        [dh, m] = NEWTON_GMRES(F,F_old,iters, M, PARAMS_LK, DIM_LK, @VERIF_FVM, h, h_old, S_old, phi_old, k_old, t);
        
        % Update estimate for current timestep's h
        h = LineSearch(DIM_LK, @VERIF_FVM, h, dh, h_old, S_old, phi_old, k_old, t, PARAMS_LK);
        h_store_LK(:,end+1)=h;
        % Update F and all other variables for this time step
        F_old = F;
        [F, S, phi, k,PUMPS, EVAPOT] = VERIF_FVM(DIM_LK, h, h_old, S_old, phi_old, k_old, t, PARAMS_LK);
        err = norm(F, 2);
        iters = iters + 1;
        Pumps_LK(:,end+1)=PUMPS;
        Evapot_LK(:,end+1)=EVAPOT;
        
        
        % If haven't converged but has been too many iterations, halve time
        % step
        if iters == PARAMS_LK.max_iters - 1 || err > 1e12
            iters = 0;
            t = t - PARAMS_LK.dt;
            PARAMS_LK.dt = PARAMS_LK.dt / 2;
            t = t + PARAMS_LK.dt;
        end
    end
    
    % Output some debug info if wanted
    if PARAMS_LK.debug == true
        fprintf('Newton converged with GMRES for t=%d in %d iterations. Err:%d dt:%d\n', ...
            t, iters, err, PARAMS_LK.dt);
    end
    
    T_LK(end + 1) = t;
    phi_avg(end + 1) = PHI_AVG(DIM_LK, phi);
    phi_true(end + 1) = phi_0;%PHI_TRUE(DIM_LK, PARAMS_LK, t, phi_0);
    
    % We have now converged, so update variables
    h_old = h;
    S_old = S;
    phi_old = phi;
    k_old = k;

    % Recalculate base error
    F = VERIF_FVM(DIM_LK, h, h, S, phi, k, t, PARAMS_LK);
    err = norm(F, 2);
    err_old = err;
    
    % If adaptive time stepping and converged quickly, increase time step
    if iters <= PARAMS_LK.gmres_max
        PARAMS_LK.dt = min(PARAMS_LK.dt * PARAMS_LK.adaptive_timestep, PARAMS_LK.max_dt);
    end
    
    if (PARAMS_LK.realtime_plot == true)
        VISUALISE(DIM_LK, h, phi, T_LK, phi_true, phi_avg);
    end
    Iteration_Time_LK(end+1)=toc;
    disp(['Iter Time ' num2str(Iteration_Time_LK(end))])
    Total_Time_LK=Total_Time_LK+Iteration_Time_LK(end);

end

disp('Steady State Reached')
disp(['Total Time ' num2str(Total_Time_LK(end))])
    
save('T_LK')
save('h_store_LK')
save('Pumps_LK')
save('Evapot_LK')
save('Iteration_Time_LK')
save('Total_Time_LK')

% CREATE_VIDEO(wcontvideo, watercontent, 20);
% CREATE_VIDEO(pheadvideo, pressurehead, 20);

