close all
clear
clc
format compact
%% Part 0 Initialisation
framenum = 1;

[PARAMS] = INIT_PARAMS;
[DIM] = GRIDCOORD;
[h_old, S_old, phi_old, k_old] = INITCOND(DIM);
% [M]=MATGEN(DIM,h_old);
% [S]=SGEN(DIM);

h = h_old;
S = S_old;
phi = phi_old;
k = k_old;

F = FVM_TEST(DIM, h, h_old, S_old, phi_old, k_old, PARAMS.dt, PARAMS);
err = norm(F, 2);
err_old = err;

iters = 0;
fevals = 0;
f_eval_total = 1;

J = jac_func_2(DIM, F, @FVM_TEST, h, h_old, S_old, phi_old, k_old, PARAMS.dt, PARAMS);


%%
videoName_wcont = 'WaterContent.avi';
videoName_phead = 'PressureHead.avi';

wcontvideo = VideoWriter(videoName_wcont);
pheadvideo = VideoWriter(videoName_phead);


if PARAMS.realtime_plot
    head_figure = figure('Name', 'Head');
    phi_figure = figure('Name', 'Water Content');
    SOL_VIS(DIM, head_figure, 'gray', ['Pressure Head (m) Time: ', num2str(0)], h);
    SOL_VIS(DIM, phi_figure, 'default', ['Water Content Time: ', num2str(0)], phi);
end

t = 0;
timesteps = 0;
T=[0];
steady_state = false;

while steady_state == false %(norm(phi-phi_old) > PARAMS.breaktol)
    t = t + PARAMS.dt;
    timesteps = timesteps + 1;
%    R=PARAMS.r_f*(1+cos((2*pi)*t/365)); For when we finally want to stop
%    using constant rainfall
    
    
    while err > PARAMS.tol_a + PARAMS.tol_r * err_old && iters < PARAMS.max_iters
        if mod(iters, PARAMS.jacobian_update) == 0
            J_old=J;
            J = jac_func_2(DIM, F, @FVM_TEST, h, h_old, S_old, phi_old, k_old, t, PARAMS);
            fevals = fevals + DIM.n * DIM.m;
        end
    
        dh = J\(-F); %This line is now in Jsolv

        
        h = LineSearch(DIM, @FVM_TEST, h, dh, h_old, S_old, phi_old, k_old, t, PARAMS);

        [F, S, phi, k] = FVM_TEST(DIM, h, h_old, S_old, phi_old, k_old, t, PARAMS);
        err = norm(F, 2);
        iters = iters + 1;
        fevals = fevals + 1;
        
        if PARAMS.debug == true
            fprintf('t:%d iters:%d err:%d fevals:%d dank:%d max_h:%d\n', t, iters, err, fevals, timesteps, max(h));
        end
        
        if iters == PARAMS.max_iters - 1 || err > 1e12
            iters = 0;
            t = t - PARAMS.dt;
            PARAMS.dt = PARAMS.dt / 2;
            t = t + PARAMS.dt;
        end
        
    end
    
    % We have now converged
    h_old = h;
    S_old = S;
    phi_old = phi;
    k_old = k;
    
    F = FVM_TEST(DIM, h, h, S, phi, k, t, PARAMS);
    err = norm(F, 2);
    err_old = err; % The error for the new timestep to compare with the next Newton iterations
    
    f_eval_total = f_eval_total + fevals + 1;
    fevals = 0;
    
    if iters <= 3
        PARAMS.dt = PARAMS.dt * 1.5;
    end
    
    iters = 0;

    if PARAMS.realtime_plot == true
        SOL_VIS(DIM, head_figure, 'gray', ['Pressure Head (m) Time: ', num2str(t)], h);
        pressurehead(framenum) = getframe(gcf);
        SOL_VIS(DIM, phi_figure, 'default', ['Water Content Time: ', num2str(t)], phi);
        watercontent(framenum) = getframe(gcf);
        framenum = framenum +1;
    end
    T(end+1)=t;
    
    if t > 100000
        steady_state = true;
    end
    
end

% CREATE_VIDEO(wcontvideo, watercontent, 20);
% CREATE_VIDEO(pheadvideo, pressurehead, 20);

PARAMS.PUMPS=1;

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







