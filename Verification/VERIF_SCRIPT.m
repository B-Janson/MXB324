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
S = S_old;
phi = phi_old;
k = k_old;

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

% videoName_wcont = 'WaterContent.avi';
% videoName_phead = 'PressureHead.avi';

% wcontvideo = VideoWriter(videoName_wcont);
% pheadvideo = VideoWriter(videoName_phead);

if PARAMS.realtime_plot
    % if realtime plotting, show the initial state of solution
    head_figure = figure('Name', 'Head');
    phi_figure = figure('Name', 'Water Content');
    SOL_VIS(DIM, head_figure, 'gray', ['Pressure Head (m) Time: ', num2str(0)], h);
    SOL_VIS(DIM, phi_figure, 'default', ['Water Content Time: ', num2str(0)], phi);
end

t = 0;
timesteps = 0;
% T=[0];
steady_state = false;

%% Part 1 Main Solver

while steady_state == false %(norm(phi-phi_old) > PARAMS.breaktol)
    t = t + PARAMS.dt;
    timesteps = timesteps + 1;
%    R=PARAMS.r_f*(1+cos((2*pi)*t/365)); For when we finally want to stop
%    using constant rainfall
    
    while err > PARAMS.tol_a + PARAMS.tol_r * err_old && iters < PARAMS.max_iters
        if mod(iters, PARAMS.jacobian_update) == 0
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
            fprintf('t:%d iters:%d err:%d fevals:%d timesteps:%d min(h):%d\n', t, iters, err, fevals, timesteps, min(h) >= 0);
        end
        
        % If haven't converged but has been too many iterations, halve time
        % step
        if iters == PARAMS.max_iters - 1 || err > 1e12
            iters = 0;
            t = t - PARAMS.dt;
            PARAMS.dt = PARAMS.dt / 2;
            t = t + PARAMS.dt;
        end
        
        % If pressure head is positive at surface, steady state reached
        if min(h) >= 0
            fprintf('steady state\n')
            steady_state = true;
            break
        end
        
    end
    
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
    if iters <= 5 && PARAMS.adaptive_timestep == true
        PARAMS.dt = PARAMS.dt * 1.2;
    end
    
    % reset iters
    iters = 0;

    if PARAMS.realtime_plot == true
        SOL_VIS(DIM, head_figure, 'gray', ['Pressure Head (m) Time: ', num2str(t)], h);
%         pressurehead(framenum) = getframe(gcf);
        SOL_VIS(DIM, phi_figure, 'default', ['Water Content Time: ', num2str(t)], phi);
%         watercontent(framenum) = getframe(gcf);
        framenum = framenum +1;
    end
%     T(end+1)=t;
    
end

% CREATE_VIDEO(wcontvideo, watercontent, 20);
% CREATE_VIDEO(pheadvideo, pressurehead, 20);

PARAMS.PUMPS=1;

