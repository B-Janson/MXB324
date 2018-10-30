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
PARAMS_BUDGET_PUMPS=PARAMS_SSL;
PARAMS_BUDGET_PUMPS.PUMPING=1;
PARAMS_BUDGET_PUMPS.endtime=50*365;
PARAMS_BUDGET_PUMPS.realtime_plot=1;
PARAMS_BUDGET_PUMPS.PUMPS=[ 450, 10, 0.213; 100, 50, 0.427];
save('PARAMS_BUDGET_PUMPS')
% Get the grid & other info about grid
load('DIM_SSL');
DIM_BUDGET_PUMPS=DIM_SSL;

%Recalculate 'initial conditions' using steady state
load('h_store_SSL');
h = h_store_SSL(:,end);
h_old=h;
h_store_BUDGET_PUMPS=h_old;
S=zeros(DIM_BUDGET_PUMPS.n*DIM_BUDGET_PUMPS.m,1);
phi=zeros(DIM_BUDGET_PUMPS.n*DIM_BUDGET_PUMPS.m,1);
k=zeros(DIM_BUDGET_PUMPS.n*DIM_BUDGET_PUMPS.m,1);
for i = 1:DIM_BUDGET_PUMPS.n*DIM_BUDGET_PUMPS.m
    S(i) = SATURATION(DIM_BUDGET_PUMPS, h, i);
    phi(i) = WATER_CONTENT(DIM_BUDGET_PUMPS, h, S, i);
    k(i) = PERM(DIM_BUDGET_PUMPS, h, S, i);
end

% Initial calculation of F at t = dt
F = VERIF_FVM(DIM_BUDGET_PUMPS, h, h_old, S_old, phi_old, k_old, PARAMS_BUDGET_PUMPS.dt, PARAMS_BUDGET_PUMPS);
err = norm(F, 2);
err_old = err;
PARAMS_BUDGET_PUMPS.F0=norm(F,inf);
framenum = 1;
F_old=F;
% Get the jacobian
J = JAC_FUNC(DIM_BUDGET_PUMPS, F, @VERIF_FVM, h, h_old, S_old, phi_old, k_old, PARAMS_BUDGET_PUMPS.dt, PARAMS_BUDGET_PUMPS);
M = J;

h = h_old;
S = S_old;
phi = phi_old;
k = k_old;

T_BUDGET_PUMPS = 0;
phi_avg = PHI_AVG(DIM_BUDGET_PUMPS, phi);
phi_0 = phi_avg;
phi_true = phi_0;

% videoName_wcont = 'WaterContent.avi';
% videoName_phead = 'PressureHead.avi';

% wcontvideo = VideoWriter(videoName_wcont);
% pheadvideo = VideoWriter(videoName_phead);

if PARAMS_BUDGET_PUMPS.realtime_plot
    % if realtime plotting, show the initial state of solution
    figure('Position', [100 160 850 500]);
    VISUALISE(DIM_BUDGET_PUMPS, h, phi, T_BUDGET_PUMPS, phi_true, phi_avg);
end

t = 0;
timesteps = 0;
steady_state = false;
m = PARAMS_BUDGET_PUMPS.gmres_max;
Pumps_BUDGET_PUMPS=zeros(DIM_BUDGET_PUMPS.n*DIM_BUDGET_PUMPS.m,1);
Evapot_BUDGET_PUMPS=zeros(DIM_BUDGET_PUMPS.n*DIM_BUDGET_PUMPS.m,1);

change_in_psi = 0;
change_in_river = 0;
change_in_rain = 0;
change_in_pumps = 0;
change_in_evapo = 0;
%% Part 1 Main Solver
Overhead_Time_BUDGET_PUMPS=toc
Iteration_Time_BUDGET_PUMPS=0;
Total_Time_BUDGET_PUMPS=Overhead_Time_BUDGET_PUMPS;
% Time step iteration
while t < PARAMS_BUDGET_PUMPS.endtime
    tic
    t = t + PARAMS_BUDGET_PUMPS.dt;
    timesteps = timesteps + 1;
    iters = 0;
    F_0_norm = norm(F, Inf);
    
    % Newton step iteration
    while err > PARAMS_BUDGET_PUMPS.tol_a + PARAMS_BUDGET_PUMPS.tol_r * err_old && iters < PARAMS_BUDGET_PUMPS.max_iters

        if m > PARAMS_BUDGET_PUMPS.gmres_max
            fprintf('Recalculating Preconditioner\n');
            M=J;
            %[M,~] = lu(J); 
            %[M,~] = ilu(J); 
            %M=diag(J); % Invariant subspace
            J = JAC_FUNC(DIM_BUDGET_PUMPS, F, @VERIF_FVM, h, h_old, S_old, phi_old, k_old, t, PARAMS_BUDGET_PUMPS);

        end
        
        % Get the del h
        [dh, m] = NEWTON_GMRES(F,F_old,iters, M, PARAMS_BUDGET_PUMPS, DIM_BUDGET_PUMPS, @VERIF_FVM, h, h_old, S_old, phi_old, k_old, t);
        
        % Update estimate for current timestep's h
        h = LineSearch(DIM_BUDGET_PUMPS, @VERIF_FVM, h, dh, h_old, S_old, phi_old, k_old, t, PARAMS_BUDGET_PUMPS);
        h_store_BUDGET_PUMPS(:,end+1)=h;
        % Update F and all other variables for this time step
        F_old = F;
        [F, S, phi, k,PUMPS, EVAPOT, river, rainfall] = VERIF_FVM(DIM_BUDGET_PUMPS, h, h_old, S_old, phi_old, k_old, t, PARAMS_BUDGET_PUMPS);
        err = norm(F, 2);
        iters = iters + 1;
        Pumps_BUDGET_PUMPS(:,end+1)=PUMPS;
        Evapot_BUDGET_PUMPS(:,end+1)=EVAPOT;
        
        
        % If haven't converged but has been too many iterations, halve time
        % step
        if iters == PARAMS_BUDGET_PUMPS.max_iters - 1 || err > 1e12
            iters = 0;
            t = t - PARAMS_BUDGET_PUMPS.dt;
            PARAMS_BUDGET_PUMPS.dt = PARAMS_BUDGET_PUMPS.dt / 2;
            t = t + PARAMS_BUDGET_PUMPS.dt;
        end
    end
    
    % Output some debug info if wanted
    if PARAMS_BUDGET_PUMPS.debug == true
        fprintf('Newton converged with GMRES for t=%d in %d iterations. Err:%d dt:%d\n', ...
            t, iters, err, PARAMS_BUDGET_PUMPS.dt);
    end
    
    T_BUDGET_PUMPS(end + 1) = t;
    phi_avg(end + 1) = PHI_AVG(DIM_BUDGET_PUMPS, phi);
    phi_true(end + 1) = phi_0;%PHI_TRUE(DIM_BUDGET_PUMPS, PARAMS_BUDGET_PUMPS, t, phi_0);
    
    change_in_psi(end + 1) = phi_avg(end) - phi_avg(end - 1);
    change_in_river(end + 1) = -river;
    change_in_rain(end + 1) = rainfall;
    change_in_pumps(end + 1) = sum(PUMPS) / (100);
    change_in_evapo(end + 1) = sum(EVAPOT) / (350 * 2 + 150 * 4);
    
    fprintf('Change in water content: %d ... Rainfall: %d ... River discharge: %s ... Pumps:%d ... Evaporation:%d\n', ...
        change_in_psi(end), rainfall, -river, change_in_pumps(end), change_in_evapo(end));
    
    fprintf('LHS:%d \t RHS%d\n', ...
        change_in_psi(end), change_in_rain(end) + change_in_river(end) + change_in_pumps(end) + change_in_evapo(end));
    
    % We have now converged, so update variables
    h_old = h;
    S_old = S;
    phi_old = phi;
    k_old = k;

    % Recalculate base error
    F = VERIF_FVM(DIM_BUDGET_PUMPS, h, h, S, phi, k, t, PARAMS_BUDGET_PUMPS);
    err = norm(F, 2);
    err_old = err;
    
    % If adaptive time stepping and converged quickly, increase time step
    if iters <= PARAMS_BUDGET_PUMPS.gmres_max
        PARAMS_BUDGET_PUMPS.dt = min(PARAMS_BUDGET_PUMPS.dt * PARAMS_BUDGET_PUMPS.adaptive_timestep, PARAMS_BUDGET_PUMPS.max_dt);
    end
    
    if (PARAMS_BUDGET_PUMPS.realtime_plot == true)
        VISUALISE(DIM_BUDGET_PUMPS, h, phi, T_BUDGET_PUMPS, phi_true, phi_avg);
    end
    Iteration_Time_BUDGET_PUMPS(end+1)=toc;
    disp(['Iter Time ' num2str(Iteration_Time_BUDGET_PUMPS(end))])
    Total_Time_BUDGET_PUMPS=Total_Time_BUDGET_PUMPS+Iteration_Time_BUDGET_PUMPS(end);

end

disp('Steady State Reached')
disp(['Total Time ' num2str(Total_Time_BUDGET_PUMPS(end))])
    
save('T_BUDGET_PUMPS')
save('h_store_BUDGET_PUMPS')
save('Pumps_BUDGET_PUMPS')
save('Evapot_BUDGET_PUMPS')
save('Iteration_Time_BUDGET_PUMPS')
save('Total_Time_BUDGET_PUMPS')

% CREATE_VIDEO(wcontvideo, watercontent, 20);
% CREATE_VIDEO(pheadvideo, pressurehead, 20);

