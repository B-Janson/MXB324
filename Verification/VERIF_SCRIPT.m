close all
clear
clc
format compact
%% Part 0 Initialisation
% Get the parameters to solve for
[PARAMS] = VERIF_PARAMS;
% Get the grid & other info about grid
[DIM] = VERIF_GRID_COORD;
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
J = JAC_FUNC(DIM, F, @VERIF_FVM, h, h_old, S_old, phi_old, k_old, PARAMS.dt, PARAMS, 'full');
% total_bandwidth = bandwidth(J);
% [DIM, h_old, S_old, phi_old, k_old, F] = REORDER_NODES(DIM, J, h_old, S_old, phi_old, k_old, F);

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

sat_col = zeros(30, 3);

for i = 1:20
    sat_col(i, 1) = max(1.0 - 1.0 * (i-1) / 10, 0);
    sat_col(i, 2) = max(1.0 - 1.0 * (i-1) / 20, 0);
    sat_col(i, 3) = 1.0;
end

for i = 21:30
    sat_col(i, 3) = 1.0 - 0.5 * (i-21) / 10;
end

if PARAMS.realtime_plot
    figure
    hold on
    plot(0, phi_avg, 'b')
    plot(0, phi_true, 'r')
    hold off
    drawnow
    % if realtime plotting, show the initial state of solution
%     head_figure = figure('Name', 'Head');
%     phi_figure = figure('Name', 'Water Content');
%     sat_figure = figure('Name', 'Saturation');
%     SOL_VIS(DIM, head_figure, 'gray', ['Pressure Head (m) Time: ', num2str(0)], h);
%     SOL_VIS(DIM, phi_figure, sat_col, ['Water Content Time: ', num2str(0)], phi);
%     SOL_VIS(DIM, sat_figure, sat_col, ['Saturation Time: ', num2str(0)], S);
end

t = 0;
timesteps = 0;
% T=[0];
steady_state = false;

%% Part 1 Main Solver
tic;
while steady_state == false && t < PARAMS.endtime %(norm(phi-phi_old) > PARAMS.breaktol)
    t = t + PARAMS.dt;
    timesteps = timesteps + 1;
%    R=PARAMS.r_f*(1+cos((2*pi)*t/365)); For when we finally want to stop
%    using constant rainfall
    
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
            PARAMS.dt = PARAMS.dt / 2;
            t = t + PARAMS.dt;
            if PARAMS.dt < 5
                PARAMS.method = 'full';
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
        hold on
        plot(T, phi_avg, 'b')
        plot(T, phi_true, 'r')
        hold off
        drawnow
%         SOL_VIS(DIM, head_figure, 'gray', ['Pressure Head (m) Time: ', num2str(t)], h);
%         pressurehead(framenum) = getframe(gcf);
%         SOL_VIS(DIM, phi_figure, sat_col, ['Water Content Time: ', num2str(t)], phi);
%         watercontent(framenum) = getframe(gcf);
%         SOL_VIS(DIM, sat_figure, sat_col, ['Saturation Time: ', num2str(t)], S);
%         framenum = framenum +1;
    end
%     T(end+1)=t;
    
end

toc

% CREATE_VIDEO(wcontvideo, watercontent, 20);
% CREATE_VIDEO(pheadvideo, pressurehead, 20);

PARAMS.PUMPS=1;

