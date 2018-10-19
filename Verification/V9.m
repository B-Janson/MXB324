function f = V9(DIM, i, h, h_old, phi, phi_old, k, k_old, t, PARAMS)
% Sandstone/Rain

n = DIM.n;

point = (DIM.r == i);
west = (DIM.r == i-1);
south = (DIM.r == i-n);

DELTA = DIM.DELTA(point, :);
dx = DELTA(1);
dz = DELTA(3);
ST = DIM.ST(point, 3);
K_xx = DIM.K_xx(ST);
K_zz = DIM.K_zz(ST);
cell_volume = DIM.VOL(point, 5);

dt = PARAMS.dt;
theta = PARAMS.theta;
r_f = PARAMS.r_f(PARAMS.r_t);
active_flow = 1;

if h_old(point) >= 0
   active_flow = 0; 
end

k_w = (k(point) + k(west)) / 2;
k_s = (k(point) + k(south)) / 2;

k_w_old = (k_old(point) + k_old(west)) / 2;
k_s_old = (k_old(point) + k_old(south)) / 2;

gamma_1 = active_flow * r_f * dx / 2;
gamma_2 = k_w * K_xx * dz / 2 * ((h(west) - h(point))/dx);
gamma_3 = k_s * K_zz * dx / 2 * ((h(south) - h(point))/dz - 1);

gamma_1_old = active_flow * r_f * dx / 2;
gamma_2_old = k_w_old * K_xx * dz / 2 * ((h_old(west) - h_old(point))/dx);
gamma_3_old = k_s_old * K_zz * dx / 2 * ((h_old(south) - h_old(point))/dz - 1);

f = phi(point) - phi_old(point) - dt/cell_volume * (theta * ...
            (gamma_1 + gamma_2 + gamma_3) ...
            + (1 - theta) * ( ...
            gamma_1_old + gamma_2_old + gamma_3_old));
        
end
