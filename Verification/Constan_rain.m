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
PARAMS_CONS=PARAMS_SSL;
PARAMS_CONS.PUMPING=1;
PARAMS_CONS.endtime=10*365;
PARAMS_CONS.realtime_plot=1;
save('PARAMS_CONS')
% Get the grid & other info about grid
load('DIM_SSL');
DIM_CONS=DIM_SSL;

%Recalculate 'initial conditions' using steady state
load('h_store_SSL');
h = h_store_SSL(:,end);
h_old=h;
h_store_CONS=h_old;
S=zeros(DIM_CONS.n*DIM_CONS.m,1);
phi=zeros(DIM_CONS.n*DIM_CONS.m,1);
k=zeros(DIM_CONS.n*DIM_CONS.m,1);
for i = 1:DIM_CONS.n*DIM_CONS.m
    S(i) = SATURATION(DIM_CONS, h, i);
    phi(i) = WATER_CONTENT(DIM_CONS, h, S, i);
    k(i) = PERM(DIM_CONS, h, S, i);
end
% Initial calculation of F at t = dt
F = VERIF_FVM(DIM_CONS, h, h_old, S_old, phi_old, k_old, PARAMS_CONS.dt, PARAMS_CONS);
err = norm(F, 2);
err_old = err;
PARAMS_CONS.F0=norm(F,inf);
framenum = 1;
F_old=F;
% Get the jacobian
J = JAC_FUNC(DIM_CONS, F, @VERIF_FVM, h, h_old, S_old, phi_old, k_old, PARAMS_CONS.dt, PARAMS_CONS);
M = J;

h = h_old;
S = S_old;
phi = phi_old;
k = k_old;

T_CONS = 0;
phi_avg = PHI_AVG(DIM_CONS, phi);
phi_0 = phi_avg;
phi_true = PHI_TRUE(DIM_CONS, PARAMS_CONS, 0, phi_0);

% videoName_wcont = 'WaterContent.avi';
% videoName_phead = 'PressureHead.avi';

% wcontvideo = VideoWriter(videoName_wcont);
% pheadvideo = VideoWriter(videoName_phead);

if PARAMS_CONS.realtime_plot
    % if realtime plotting, show the initial state of solution
    figure('Position', [100 160 850 500]);
    VISUALISE(DIM_CONS, h, phi, T_CONS, phi_true, phi_avg);
end

t = 0;
timesteps = 0;
steady_state = false;
m = PARAMS_CONS.gmres_max;
Pumps_CONS=zeros(DIM_CONS.n*DIM_CONS.m,1);
Evapot_CONS=zeros(DIM_CONS.n*DIM_CONS.m,1);
%% Part 1 Main Solver
Overhead_Time_CONS=toc
Iteration_Time_CONS=0;
Total_Time_CONS=Overhead_Time_CONS;
% Time step iteration
while t < PARAMS_CONS.endtime
    tic
    t = t + PARAMS_CONS.dt;
    timesteps = timesteps + 1;
    iters = 0;
    F_0_norm = norm(F, Inf);
    
    % Newton step iteration
    while err > PARAMS_CONS.tol_a + PARAMS_CONS.tol_r * err_old && iters < PARAMS_CONS.max_iters

        if m > PARAMS_CONS.gmres_max
            fprintf('Recalculating Preconditioner\n');
            M=J;
            %[M,~] = lu(J); 
            %[M,~] = ilu(J); 
            %M=diag(J); % Invariant subspace
            J = JAC_FUNC(DIM_CONS, F, @VERIF_FVM, h, h_old, S_old, phi_old, k_old, t, PARAMS_CONS);

        end
        
        % Get the del h
        [dh, m] = NEWTON_GMRES(F,F_old,iters, M, PARAMS_CONS, DIM_CONS, @VERIF_FVM, h, h_old, S_old, phi_old, k_old, t);
        
        % Update estimate for current timestep's h
        h = LineSearch(DIM_CONS, @VERIF_FVM, h, dh, h_old, S_old, phi_old, k_old, t, PARAMS_CONS);
        h_store_CONS(:,end+1)=h;
        % Update F and all other variables for this time step
        F_old = F;
        [F, S, phi, k,PUMPS, EVAPOT] = VERIF_FVM(DIM_CONS, h, h_old, S_old, phi_old, k_old, t, PARAMS_CONS);
        err = norm(F, 2);
        iters = iters + 1;
        Pumps_CONS(:,end+1)=PUMPS;
        Evapot_CONS(:,end+1)=EVAPOT;
        
        
        % If haven't converged but has been too many iterations, halve time
        % step
        if iters == PARAMS_CONS.max_iters - 1 || err > 1e12
            iters = 0;
            t = t - PARAMS_CONS.dt;
            PARAMS_CONS.dt = PARAMS_CONS.dt / 2;
            t = t + PARAMS_CONS.dt;
        end
    end
    
    % Output some debug info if wanted
    if PARAMS_CONS.debug == true
        fprintf('Newton converged with GMRES for t=%d in %d iterations. Err:%d dt:%d\n', ...
            t, iters, err, PARAMS_CONS.dt);
    end
    
    T_CONS(end + 1) = t;
    phi_avg(end + 1) = PHI_AVG(DIM_CONS, phi);
    phi_true(end + 1) = PHI_TRUE(DIM_CONS, PARAMS_CONS, t, phi_0);
    
    % We have now converged, so update variables
    h_old = h;
    S_old = S;
    phi_old = phi;
    k_old = k;

    % Recalculate base error
    F = VERIF_FVM(DIM_CONS, h, h, S, phi, k, t, PARAMS_CONS);
    err = norm(F, 2);
    err_old = err;
    
    % If adaptive time stepping and converged quickly, increase time step
    if iters <= PARAMS_CONS.gmres_max
        PARAMS_CONS.dt = min(PARAMS_CONS.dt * PARAMS_CONS.adaptive_timestep, PARAMS_CONS.max_dt);
    end
    
    if (PARAMS_CONS.realtime_plot == true)
        VISUALISE(DIM_CONS, h, phi, T_CONS, phi_true, phi_avg);
    end
    Iteration_Time_CONS(end+1)=toc;
    disp(['Iter Time ' num2str(Iteration_Time_CONS(end))])
    Total_Time_CONS=Total_Time_CONS+Iteration_Time_CONS(end);

end

disp('Steady State Reached')
disp(['Total Time ' num2str(Total_Time_CONS(end))])
    
save('T_CONS')
save('h_store_CONS')
save('Pumps_CONS')
save('Evapot_CONS')
save('Iteration_Time_CONS')
save('Total_Time_CONS')

% CREATE_VIDEO(wcontvideo, watercontent, 20);
% CREATE_VIDEO(pheadvideo, pressurehead, 20);

