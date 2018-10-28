%% Constant Rainfall,higher conductivity With Pumps on
%Run from the steady state for 10 years to invesigate pumping


close all
clear
clc
format compact
%% Part 0 Initialisation
% Get the parameters to solve for
tic
 load('PARAMS_SSL')
PARAMS_HK=PARAMS_SSL;
PARAMS_HK.PUMPING=1;
PARAMS_HK.endtime=25*365;
PARAMS_HK.realtime_plot=1;
PARAMS_HK.K_r=0.6;
PARAMS_HK.PUMPS=[ 450, 10, 0.213; 100, 50, 0.427];
save('PARAMS_HK')
% Get the grid & other info about grid
load('DIM_SSL');
DIM_HK=DIM_SSL;

%Recalculate 'initial conditions' using steady state
load('h_store_SSL');
h = h_store_SSL(:,end);
h_old=h;
h_store_HK=h_old;
S=zeros(DIM_HK.n*DIM_HK.m,1);
phi=zeros(DIM_HK.n*DIM_HK.m,1);
k=zeros(DIM_HK.n*DIM_HK.m,1);
for i = 1:DIM_HK.n*DIM_HK.m
    S(i) = SATURATION(DIM_HK, h, i);
    phi(i) = WATER_CONTENT(DIM_HK, h, S, i);
    k(i) = PERM(DIM_HK, h, S, i);
end

% Initial calculation of F at t = dt
F = VERIF_FVM(DIM_HK, h, h_old, S_old, phi_old, k_old, PARAMS_HK.dt, PARAMS_HK);
err = norm(F, 2);
err_old = err;
PARAMS_HK.F0=norm(F,inf);
framenum = 1;
F_old=F;
% Get the jacobian
J = JAC_FUNC(DIM_HK, F, @VERIF_FVM, h, h_old, S_old, phi_old, k_old, PARAMS_HK.dt, PARAMS_HK);
M = J;

h = h_old;
S = S_old;
phi = phi_old;
k = k_old;

T_HK = 0;
phi_avg = PHI_AVG(DIM_HK, phi);
phi_0 = phi_avg;
phi_true = phi_0;

% videoName_wcont = 'WaterContent.avi';
% videoName_phead = 'PressureHead.avi';

% wcontvideo = VideoWriter(videoName_wcont);
% pheadvideo = VideoWriter(videoName_phead);

if PARAMS_HK.realtime_plot
    % if realtime plotting, show the initial state of solution
    figure('Position', [100 160 850 500]);
    VISUALISE(DIM_HK, h, phi, T_HK, phi_true, phi_avg);
end

t = 0;
timesteps = 0;
steady_state = false;
m = PARAMS_HK.gmres_max;
Pumps_HK=zeros(DIM_HK.n*DIM_HK.m,1);
Evapot_HK=zeros(DIM_HK.n*DIM_HK.m,1);
%% Part 1 Main Solver
Overhead_Time_HK=toc
Iteration_Time_HK=0;
Total_Time_HK=Overhead_Time_HK;
% Time step iteration
while t < PARAMS_HK.endtime
    tic
    t = t + PARAMS_HK.dt;
    timesteps = timesteps + 1;
    iters = 0;
    F_0_norm = norm(F, Inf);
    
    % Newton step iteration
    while err > PARAMS_HK.tol_a + PARAMS_HK.tol_r * err_old && iters < PARAMS_HK.max_iters

        if m > PARAMS_HK.gmres_max
            fprintf('Recalculating Preconditioner\n');
            M=J;
            %[M,~] = lu(J); 
            %[M,~] = ilu(J); 
            %M=diag(J); % Invariant subspace
            J = JAC_FUNC(DIM_HK, F, @VERIF_FVM, h, h_old, S_old, phi_old, k_old, t, PARAMS_HK);

        end
        
        % Get the del h
        [dh, m] = NEWTON_GMRES(F,F_old,iters, M, PARAMS_HK, DIM_HK, @VERIF_FVM, h, h_old, S_old, phi_old, k_old, t);
        
        % Update estimate for current timestep's h
        h = LineSearch(DIM_HK, @VERIF_FVM, h, dh, h_old, S_old, phi_old, k_old, t, PARAMS_HK);
        h_store_HK(:,end+1)=h;
        % Update F and all other variables for this time step
        F_old = F;
        [F, S, phi, k,PUMPS, EVAPOT] = VERIF_FVM(DIM_HK, h, h_old, S_old, phi_old, k_old, t, PARAMS_HK);
        err = norm(F, 2);
        iters = iters + 1;
        Pumps_HK(:,end+1)=PUMPS;
        Evapot_HK(:,end+1)=EVAPOT;
        
        
        % If haven't converged but has been too many iterations, halve time
        % step
        if iters == PARAMS_HK.max_iters - 1 || err > 1e12
            iters = 0;
            t = t - PARAMS_HK.dt;
            PARAMS_HK.dt = PARAMS_HK.dt / 2;
            t = t + PARAMS_HK.dt;
        end
    end
    
    % Output some debug info if wanted
    if PARAMS_HK.debug == true
        fprintf('Newton converged with GMRES for t=%d in %d iterations. Err:%d dt:%d\n', ...
            t, iters, err, PARAMS_HK.dt);
    end
    
    T_HK(end + 1) = t;
    phi_avg(end + 1) = PHI_AVG(DIM_HK, phi);
    phi_true(end + 1) = phi_0;%PHI_TRUE(DIM_HK, PARAMS_HK, t, phi_0);
    
    % We have now converged, so update variables
    h_old = h;
    S_old = S;
    phi_old = phi;
    k_old = k;

    % Recalculate base error
    F = VERIF_FVM(DIM_HK, h, h, S, phi, k, t, PARAMS_HK);
    err = norm(F, 2);
    err_old = err;
    
    % If adaptive time stepping and converged quickly, increase time step
    if iters <= PARAMS_HK.gmres_max
        PARAMS_HK.dt = min(PARAMS_HK.dt * PARAMS_HK.adaptive_timestep, PARAMS_HK.max_dt);
    end
    
    if (PARAMS_HK.realtime_plot == true)
        VISUALISE(DIM_HK, h, phi, T_HK, phi_true, phi_avg);
    end
    Iteration_Time_HK(end+1)=toc;
    disp(['Iter Time ' num2str(Iteration_Time_HK(end))])
    Total_Time_HK=Total_Time_HK+Iteration_Time_HK(end);

end

disp('Steady State Reached')
disp(['Total Time ' num2str(Total_Time_HK(end))])
    
save('T_HK')
save('h_store_HK')
save('Pumps_HK')
save('Evapot_HK')
save('Iteration_Time_HK')
save('Total_Time_HK')

% CREATE_VIDEO(wcontvideo, watercontent, 20);
% CREATE_VIDEO(pheadvideo, pressurehead, 20);

