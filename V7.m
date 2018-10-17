function f = V7(DIM, i, h, h_old, phi, phi_old, k, k_old, t, PARAMS)
% Top Left

n = DIM.n;

point = (DIM.r == i);
east = (DIM.r == i+1);
south = (DIM.r == i-n);

DELTA = DIM.DELTA(point, :);
dx = DELTA(2);
dz = DELTA(3);
K_xx = DIM.K_xx(point, 4);
K_zz = DIM.K_zz(point, 4);
cell_volume = DIM.VOL(point, 5);

dt = PARAMS.dt;
theta = PARAMS.theta;
r_f = PARAMS.r_f(PARAMS.r_t);
active_flow = 1;

if h_old(point) >= 0
   active_flow = 0; 
end

k_e = (k(point) + k(east)) / 2;
k_s = (k(point) + k(south)) / 2;

k_e_old = (k_old(point) + k_old(east)) / 2;
k_s_old = (k_old(point) + k_old(south)) / 2;

gamma_1 = k_e * K_xx * dz / 2 * ((h(east) - h(point))/dx);
gamma_2 = active_flow * r_f * dx / 2;
gamma_3 = k_s * K_zz * dx / 2 * ((h(south) - h(point))/dz - 1);

gamma_1_old = k_e_old * K_xx * dz / 2 * ((h_old(east) - h_old(point))/dx);
gamma_2_old = active_flow * r_f * dx / 2;
gamma_3_old = k_s_old * K_zz * dx / 2 * ((h_old(south) - h_old(point))/dz - 1);

f = phi(point) - phi_old(point) - dt/cell_volume * ... 
            (theta * (gamma_1 + gamma_2 + gamma_3) + ...
            (1 - theta) * (gamma_1_old + gamma_2_old + gamma_3_old));
        
end
