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



jacobian = jac_func(DIM, F, @FVM_TEST, h, h_old, S_old, phi_old, k_old, PARAMS.dt, PARAMS);
total_bw = 2*bandwidth(jacobian)+1;

head_figure = figure('Name', 'Head');
phi_figure = figure('Name', 'Water Content');
SOL_VIS(DIM, head_figure, 'gray', ['Pressure Head (m) Time: ', num2str(0)], h);
SOL_VIS(DIM, phi_figure, 'default', ['Water Content Time: ', num2str(0)], phi);

t = 0;
dank = 0;

while t < PARAMS.endtime
    t = t + PARAMS.dt;
    
%     for i = 1:n*m
%         S(i) = SATURATION(DIM, h, i);
%         phi(i) = WATER_CONTENT(DIM, h, S, i);
%         k(i) = PERM(DIM, h, S, i);
%     end
    
    dank = dank + 1;
    while err > PARAMS.tol_a + PARAMS.tol_r * err_old && iters < PARAMS.max_iters
        if mod(iters, PARAMS.jacobian_update) == 0
            jacobian = jac_func(DIM, F, @FVM_TEST, h, h_old, S_old, phi_old, k_old, t, PARAMS);
            fevals = fevals + DIM.n * DIM.m;
        end
    
        del_h = jacobian\(-F);
%         h = h + del_h;
        h = lineSearch(DIM, @FVM_TEST, h, del_h, h_old, S_old, phi_old, k_old, t, PARAMS);

        F = FVM_TEST(DIM, h, h_old, S_old, phi_old, k_old, t, PARAMS);
        err = norm(F, 2);
        iters = iters + 1;
        fevals = fevals + 1;
        
        if PARAMS.debug == true
            fprintf('t:%d iters:%d err:%d fevals:%d dank:%d max_h:%d\n', t, iters, err, fevals, dank, max(h));
        end
        
        if iters == PARAMS.max_iters - 1 || err > 1e12
            iters = 0;
            t = t - PARAMS.dt;
            PARAMS.dt = PARAMS.dt / 2;
            t = t + PARAMS.dt;
        end
        
    end
    
    h_norm = norm(h_old - h, 2)
    
    h_old = h;
    
    
    for i = 1:DIM.n*DIM.m
        S(i) = SATURATION(DIM, h, i);
        phi(i) = WATER_CONTENT(DIM, h, S, i);
%         k(i) = PERM(DIM, h, S, i);
    end
    
    if iters <= 3
        PARAMS.dt = PARAMS.dt * 1.5;
    end

    if PARAMS.realtime_plot == true
        if 0 == 0
            SOL_VIS(DIM, head_figure, 'gray', ['Pressure Head (m) Time: ', num2str(t)], h);
            SOL_VIS(DIM, phi_figure, 'default', ['Water Content Time: ', num2str(t)], phi);
        end
    end
    
    
    iters = 0;
    F = FVM_TEST(DIM, h, h_old, S_old, phi_old, k_old, t, PARAMS);
    err = 1e10;
    err_old = err; % The error for the new timestep to compare with the next Newton iterations
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

%% Linesearch
function h_star = lineSearch(DIM, FVM_func, h, del_h, h_old, S_old, phi_old, k_old, t, PARAMS)


lambda = 1;
min_lambda = 10^-6;
h_star = h + lambda * del_h;

base_error = norm(FVM_func(DIM, h, h_old, S_old, phi_old, k_old, t, PARAMS));
current_error = norm(FVM_TEST(DIM, h_star, h_old, S_old, phi_old, k_old, t, PARAMS));

% fprintf('Line searching ENABLED')

while current_error > base_error && lambda > min_lambda
    lambda = lambda/2;
    h_star = h + lambda * del_h;
    current_error = norm(FVM_func(DIM, h_star, h_old, S_old, phi_old, k_old, t, PARAMS));
%     fprintf("Line search ENABLED \n")
end

if lambda < min_lambda
    fprintf("min lambda reached :)\n")
end



end


%%

% H_array = zeros(DIM.n*DIM.m, length(PARAMS.plot_times));
% 
% if PARAMS.plot_times(1) == 0
%     H_array(:, 1) = h;
%     h_count = 2;
% else
%     h_count = 1;
% end







