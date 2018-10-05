function f = V33(DIM, i, h, h_old, phi, phi_old, k, k_old, t, PARAMS)
% Sandstone/Rain

n = DIM.n;
dx = DIM.dx(1,1);
dz = DIM.dz(1,1);
K_xx = DIM.K_xx(1,1);
K_zz = DIM.K_zz(1,1);
cell_volume = DIM.cell_volume(1,1);

point = (DIM.r == i);
west = (DIM.r == i-1);
south = (DIM.r == i-n);

dt = PARAMS.dt;
theta = PARAMS.theta;
r_f = PARAMS.r_f;

k_s = 3.2e-5;

gamma_1 = -1 * r_f * dx / 2;
gamma_2 = k_s * K_zz * ((h(south) - h(point)) / dz + 1) * dx / 2;

gamma_1_old = -1 * r_f * dx / 2;
gamma_2_old = k_s * K_zz * ((h_old(south) - h_old(point)) / dz + 1) * dx / 2;

f = phi(point) - phi_old(point) + dt/cell_volume * (theta * ( ...
            gamma_1 + gamma_2) ...
            + (1 - theta) * ( ...
            gamma_1_old + gamma_2_old));
        
end
