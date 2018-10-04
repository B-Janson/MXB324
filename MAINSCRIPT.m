close all
clear
clc
%% Part 0 Initialisation
[PARAMS] = INIT_PARAMS;
[DIM] = GRIDCOORD;
[h_old, phi_old, k_old, S_old] = INITCOND(DIM);
% [M]=MATGEN(DIM,h_old);
% [S]=SGEN(DIM);

h = h_old;
phi = phi_old;
k = k_old;
S = S_old;

F = FVM_TEST(DIM, h, h_old, S_old, phi_old, k_old, PARAMS.dt, PARAMS);
err = norm(F, 2);
err_old = err;

iters = 0;
fevals = 0;
f_eval_total = 1;

H_array = zeros(DIM.n*DIM.m, length(PARAMS.plot_times));

if PARAMS.plot_times(1) == 0
    H_array(:, 1) = h;
    h_count = 2;
else
    h_count = 1;
end

jacobian = jac_func(DIM, F, @FVM_TEST, h, h_old, S_old, phi_old, k_old, 0, PARAMS);
total_bw = 2*bandwidth(jacobian)+1;

head_figure = figure('Name', 'Head');
phi_figure = figure('Name', 'Water Content');

for t = PARAMS.T
    while err > PARAMS.tol_a + PARAMS.tol_r * err_old && iters < PARAMS.max_iters
        if mod(iters, PARAMS.jacobian_update) == 0
            jacobian = jac_func(DIM, F, @FVM_TEST, h, h_old, S_old, phi_old, k_old, t, PARAMS);
            fevals = fevals + DIM.n * DIM.m;
        end
    
        del_h = jacobian\(-F);
        h = h + del_h;

        F = FVM_TEST(DIM, h, h_old, S_old, phi_old, k_old, t, PARAMS);
        err = norm(F, 2);
        iters = iters + 1;
        fevals = fevals + 1;
        
        if PARAMS.debug == true
            fprintf('t:%d iters:%d err:%d fevals:%d\n', t, iters, err, fevals);
        end
    end

    if PARAMS.realtime_plot == true
        SOL_VIS(DIM, head_figure, 'gray', ['Pressure Head (m) Time: ', num2str(t)], h);
        SOL_VIS(DIM, phi_figure, 'default', ['Water Content Time: ', num2str(t)], phi);
    end
    
    h_old = h;
    iters = 0;
    err_old = err; % The error for the new timestep to compare with the next Newton iterations
    err = 1e10;
    f_eval_total = f_eval_total + fevals + 1;
    fevals = 0;
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
