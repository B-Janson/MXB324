close all
clear
clc
format compact
%% Part 1 Initialisation
% Get the parameters to solve for
[PARAMS] = INIT_PARAMS;
% Get the grid & other info about grid
[DIM] = GRIDCOORD;
% Get initial conditions
[h_old, S_old, phi_old, k_old] = INITCOND(DIM);

h = h_old;

% Initial calculation of F at t = dt
F = FVM_SOLVE(DIM, h, h_old, S_old, phi_old, k_old, PARAMS.dt, PARAMS);
err = norm(F, 2);
err_old = err;

iters = 0;
fevals = 0;
f_eval_total = 1;
framenum = 1;

J = jac_func_2(DIM, F, @FVM_SOLVE, h, h_old, S_old, phi_old, k_old, PARAMS.dt, PARAMS, 'full');
total_bandwidth = bandwidth(J);
[DIM, h_old, S_old, phi_old, k_old, F] = REORDER_NODES(DIM, J, h_old, S_old, phi_old, k_old, F);

h = h_old;
S = S_old;
phi = phi_old;
k = k_old;

% videoName_wcont = 'WaterContent.avi';
% videoName_phead = 'PressureHead.avi';

% wcontvideo = VideoWriter(videoName_wcont);
% pheadvideo = VideoWriter(videoName_phead);

if PARAMS.realtime_plot
    head_figure = figure('Name', 'Head');
    phi_figure = figure('Name', 'Water Content');
    SOL_VIS(DIM, head_figure, 'gray', ['Pressure Head (m) Time: ', num2str(0)], h);
    SOL_VIS(DIM, phi_figure, 'default', ['Water Content Time: ', num2str(0)], phi);
end

t = 0;
timesteps = 0;
% T=[0];
steady_state = false;

%% Part 2 Main Solver (to steady state)

while steady_state == false && t < PARAMS.endtime %(norm(phi-phi_old) > PARAMS.breaktol)
    t = t + PARAMS.dt;
    timesteps = timesteps + 1;
%    R=PARAMS.r_f*(1+cos((2*pi)*t/365)); For when we finally want to stop
%    using constant rainfall
    
    while err > PARAMS.tol_a + PARAMS.tol_r * err_old && iters < PARAMS.max_iters
        rho = err / err_old;
        if mod(iters, PARAMS.jacobian_update) == 0 && err > 1e-8 || rho > PARAMS.rho_max
            J_old=J;
            J = jac_func_2(DIM, F, @FVM_SOLVE, h, h_old, S_old, phi_old, k_old, t, PARAMS);
            fevals = fevals + DIM.n * DIM.m;
        end
    
        dh = J\(-F); %This line is now in Jsolv

        h = LineSearch(DIM, @FVM_SOLVE, h, dh, h_old, S_old, phi_old, k_old, t, PARAMS);

        [F, S, phi, k] = FVM_SOLVE(DIM, h, h_old, S_old, phi_old, k_old, t, PARAMS);
        err = norm(F, 2);
        iters = iters + 1;
        fevals = fevals + 1;
        
        if PARAMS.debug == true
            fprintf('t:%d iters:%d err:%d fevals:%d timesteps:%d steady_state:%d dt:%d rho:%d method:%s norm(phi-phi_old):%d\n', ...
                t, iters, err, fevals, timesteps, min(h) >= 0, PARAMS.dt, rho, PARAMS.method, norm(phi - phi_old, 2));
        end
        
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
    
    % We have now converged
    h_old = h;
    S_old = S;
    phi_old = phi;
    k_old = k;
    
    F = FVM_SOLVE(DIM, h, h, S, phi, k, t, PARAMS);
    err = norm(F, 2);
    err_old = err; % The error for the new timestep to compare with the next Newton iterations
    
    f_eval_total = f_eval_total + fevals + 1;
    fevals = 0;
    
    if iters <= PARAMS.jacobian_update * 1.5
        PARAMS.dt = min(PARAMS.dt * PARAMS.adaptive_timestep, PARAMS.max_dt);
    end
    
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

save('steady_state');

% CREATE_VIDEO(wcontvideo, watercontent, 20);
% CREATE_VIDEO(pheadvideo, pressurehead, 20);

%% Part 3 Pumps
clear, close all
load('steady_state.mat')
PARAMS.PUMPS = 1;
% PARAMS.dt = 1;
PARAMS.endtime = 20000;

F = FVM_SOLVE(DIM, h, h, S, phi, k, t, PARAMS);
err = norm(F, 2);
err_old = err; % The error for the new timestep to compare with the next Newton iterations

while t < PARAMS.endtime %(norm(phi-phi_old) > PARAMS.breaktol)
    t = t + PARAMS.dt;
    timesteps = timesteps + 1;
%    R=PARAMS.r_f*(1+cos((2*pi)*t/365)); For when we finally want to stop
%    using constant rainfall
    
    while err > PARAMS.tol_a + PARAMS.tol_r * err_old && iters < PARAMS.max_iters
        rho = err / err_old;
        if mod(iters, PARAMS.jacobian_update) == 0 && err > 1e-8 || rho > PARAMS.rho_max
            J_old=J;
            J = jac_func_2(DIM, F, @FVM_SOLVE, h, h_old, S_old, phi_old, k_old, t, PARAMS);
            fevals = fevals + DIM.n * DIM.m;
        end
    
        dh = J\(-F); %This line is now in Jsolv

        h = LineSearch(DIM, @FVM_SOLVE, h, dh, h_old, S_old, phi_old, k_old, t, PARAMS);

        [F, S, phi, k] = FVM_SOLVE(DIM, h, h_old, S_old, phi_old, k_old, t, PARAMS);
        err = norm(F, 2);
        iters = iters + 1;
        fevals = fevals + 1;
        
        if PARAMS.debug == true
            fprintf('t:%d iters:%d err:%d fevals:%d timesteps:%d steady_state:%d dt:%d rho:%d method:%s norm(dh,2):%d\n', ...
                t, iters, err, fevals, timesteps, min(h) >= 0, PARAMS.dt, rho, PARAMS.method, norm(dh,2));
        end
        
        if iters == PARAMS.max_iters - 1 || err > 1e4
            iters = 0;
            t = t - PARAMS.dt;
            PARAMS.dt = PARAMS.dt / 2;
            t = t + PARAMS.dt;
            
            if PARAMS.dt < 5
                PARAMS.method = 'full';
            end
        end
        
    end
    
    % We have now converged
    h_old = h;
    S_old = S;
    phi_old = phi;
    k_old = k;
    
    F = FVM_SOLVE(DIM, h, h, S, phi, k, t, PARAMS);
    err = norm(F, 2);
    err_old = err; % The error for the new timestep to compare with the next Newton iterations
    
    f_eval_total = f_eval_total + fevals + 1;
    fevals = 0;
    
    if iters <= PARAMS.jacobian_update * 1.5
        PARAMS.dt = PARAMS.dt * PARAMS.adaptive_timestep;
    end
    
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

% F(DIM.r) = F;
% h(DIM.r) = h;
% phi(DIM.r) = phi;
% k(DIM.r) = k;
% S(DIM.r) = S;


% HEADVIS1(DIM, h)
% HEADVIS2(DIM, phi)

%% PART 1 Convergence?







%Storage
%U is an (m*n)x1 matrix stroing the initial condition
%
%
%for start to finish
%u(:,end+1)=last solution
%Solve system
%T(i)=timesteptaken this step
%end


%%

% H_array = zeros(DIM.n*DIM.m, length(PARAMS.plot_times));
% 
% if PARAMS.plot_times(1) == 0
%     H_array(:, 1) = h;
%     h_count = 2;
% else
%     h_count = 1;
% end







