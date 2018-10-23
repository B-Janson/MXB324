close all
clear
clc
format compact
%% Part 0 Initialisation
% Get the parameters to solve for
[PARAMS] = VERIF_PARAMS;
% Get the grid & other info about grid
% Decide on source terms

%Switch for linearity. 1=linear other=other. Like[X,Z]
LIN=[1,1];

%PUMPS. Add more lines as [x,z,fraction of rainfall]
PUMPS=[100,50,0.5;... %Town Pump
      450,10,0.25];   %Farm Pump
  
%Evapotranspiration zones. define as [L,R,Depth,fraction of rainfall]
EVAPOT=[0,350,2,0.025;... %Alluviam zone
        350,500,4,0.035]; %Sandstone zone

%Precalculate everything to do with the domain
[DIM] = VERIF_GRID_COORD(LIN,PUMPS,EVAPOT);
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
% total_bandwidth = bandwidth(J);
J2=jac_funcOG(DIM, F, @VERIF_FVM, h, h_old, S_old, phi_old, k_old, PARAMS.dt, PARAMS);
norm(J-J2)

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
    head_figure = figure('Name', 'Head');
    phi_figure = figure('Name', 'Water Content');
%     sat_figure = figure('Name', 'Saturation');
    SOL_VIS(DIM, head_figure, 'gray', ['Pressure Head (m) Time: ', num2str(0)], h);
    SOL_VIS(DIM, phi_figure, 'parula', ['Water Content Time: ', num2str(0)], phi);
%     SOL_VIS(DIM, sat_figure, sat_col, ['Saturation Time: ', num2str(0)], S);
    
    analytic_figure = figure('Name', 'Analytic');
    axis([0 PARAMS.endtime 0 0.3])
    SOL_ANALYTIC(analytic_figure, T, phi_avg, phi_true)
end

t = 0;
timesteps = 0;
steady_state = false;

%% Part 1 Main Solver
tic;
while steady_state == false && t < PARAMS.endtime
    t = t + PARAMS.dt;
    timesteps = timesteps + 1;
    
    while err > PARAMS.tol_a + PARAMS.tol_r * err_old && iters < PARAMS.max_iters
        rho = err / err_old;
        if mod(iters, PARAMS.jacobian_update) == 0 && err > 1e-8 || rho > PARAMS.rho_max
            J_old=J;
            J = JAC_FUNC(DIM, F, @VERIF_FVM, h, h_old, S_old, phi_old, k_old, t, PARAMS);
            fevals = fevals + DIM.n * DIM.m;
        end
        
        % Get the del h
        dh = J\(-F); %This line is now in Jsolv
        % Update estimate for current timestep's h
        h = LineSearch(DIM, @VERIF_FVM, h, dh, h_old, S_old, phi_old, k_old, t, PARAMS);

        % Update F and all other variables for this time step
        [F, S, phi, k] = VERIF_FVM(DIM, h, h_old, S_old, phi_old, k_old, t, PARAMS);
        err = norm(F, 2);
        iters = iters + 1;
        fevals = fevals + 1;
        disp(PARAMS.dt)
        % Output some debug info if wanted
        if PARAMS.debug == true
            fprintf('t:%d iters:%d err:%d fevals:%d timesteps:%d steady_state:%d dt:%d rho:%d method:%s\n', ...
                t, iters, err, fevals, timesteps, min(h) >= 0, PARAMS.dt, rho, PARAMS.method);
        end
        
        % If haven't converged but has been too many iterations, halve time
        % step
        if iters == PARAMS.max_iters - 1 || err > 1e12
            iters = 0;
            t = t - PARAMS.dt;
            PARAMS.dt = PARAMS.dt / 3;
            t = t + PARAMS.dt;
            if PARAMS.dt < 5
                %This should be gmres eventually
            end
        end
        
        % If pressure head is positive at surface, steady state reached
        if min(h) >= 0
            fprintf('steady state\n')
            steady_state = true;
            break
        end
        
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
    if iters <= PARAMS.jacobian_update * 1.5
        PARAMS.dt = min(PARAMS.dt * PARAMS.adaptive_timestep, PARAMS.max_dt);
    end
    
    % reset iters
    iters = 0;

    if PARAMS.realtime_plot == true
        SOL_VIS(DIM, head_figure, 'gray', ['Pressure Head (m) Time: ', num2str(t)], h);
%         pressurehead(framenum) = getframe(gcf);
        SOL_VIS(DIM, phi_figure, 'parula', ['Water Content Time: ', num2str(t)], phi);
%         watercontent(framenum) = getframe(gcf);
%         SOL_VIS(DIM, sat_figure, sat_col, ['Saturation Time: ', num2str(t)], S);
        
        SOL_ANALYTIC(analytic_figure, T, phi_avg, phi_true)
%         framenum = framenum + 1;
    end
    
end
disp('Steady State Reached')
toc

% CREATE_VIDEO(wcontvideo, watercontent, 20);
% CREATE_VIDEO(pheadvideo, pressurehead, 20);

PARAMS.PUMPS=1;

