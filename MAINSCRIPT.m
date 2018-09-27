close all
clear
clc
%% Part 0 Initialisation
T=linspace(0,24,25); %Units, maybe weeks?
dt = 0.5;
tol_a = 1e-10;
tol_r = 0;
max_iters = 100;
jacobian_update = 1;
endtime = 100;

[DIM] = GRIDCOORD;
[h_old, phi_old, k_old, S_old] = INITCOND(DIM);
% [M]=MATGEN(DIM,h_old);
%[S]=SGEN(DIM);

h = h_old;
phi = phi_old;
k = k_old;
S = S_old;
theta = 1/2;

F = FVM_TEST(DIM, h, h_old, S_old, phi_old, k_old, dt, dt, theta);
err = norm(F, 2);
err_old = err;

iters = 0;
fevals = 0;
f_eval_total = 1;

jacobian = jac_func(DIM, F, @FVM_TEST, h, h_old, S_old, phi_old, k_old, 0, dt, theta);
total_bw = 2*bandwidth(jacobian)+1;

for t = dt:dt:endtime
    while err > tol_a + tol_r * err_old && iters < max_iters
        if mod(iters, jacobian_update) == 0
            jacobian = jac_func(DIM, F, @FVM_TEST, h, h_old, S_old, phi_old, k_old, t, dt, theta);
            fevals = fevals + DIM.n * DIM.m;
        end
    
        del_h = jacobian\(-F);
        h = h + del_h;

        F = FVM_TEST(DIM, h, h_old, S_old, phi_old, k_old, t, dt, theta);
        err = norm(F, 2);
        iters = iters + 1;
        fevals = fevals + 1;
        disp(t)
    end
    
    h_old = h;
    iters = 0;
    F = FVM_TEST(DIM, h, h_old, S_old, phi_old, k_old, dt, dt, theta);
    err = norm(F,2); % Get the new error
    err_old = err; % The error for the new timestep to compare with the next Newton iterations
    f_eval_total = f_eval_total + fevals + 1;
    fevals = 0;
end

% F(DIM.r) = F;
% h(DIM.r) = h;
% phi(DIM.r) = phi;
% k(DIM.r) = k;
% S(DIM.r) = S;


HEADVIS1(DIM, h)
HEADVIS2(DIM, phi)

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
