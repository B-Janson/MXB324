function f = V1(DIM, h, h_old, phi, phi_old, k, k_old, PARAMS)
% Bottom Left Corner

n = DIM.n;

point = (DIM.r == 1);
east = (DIM.r == 2);
north = (DIM.r == n+1);

DELTA = DIM.DELTA(point, :);
dx = DELTA(2);
dz = DELTA(4);
K_xx = DIM.K_xx(point, 1);
K_zz = DIM.K_zz(point, 1);
cell_volume = DIM.VOL(point, 5);

dt = PARAMS.dt;
theta = PARAMS.theta;

k_e = (k(point) + k(east)) / 2;
k_n = (k(point) + k(north)) / 2;

k_e_old = (k_old(point) + k_old(east)) / 2;
k_n_old = (k_old(point) + k_old(north)) / 2;

gamma_1 = -k_e * K_xx * dz / 2 * ((h(east) - h(point))/dx);
gamma_2 = -k_n * K_zz * dx / 2 * (1 + (h(north) - h(point))/dz);

gamma_1_old = -k_e_old * K_xx * dz / 2 * ((h_old(east) - h_old(point))/dx);
gamma_2_old = -k_n_old * K_zz * dx / 2 * (1 + (h_old(north) - h_old(point))/dz);

f = phi(point) - phi_old(point) + dt/cell_volume * ( ...
              theta * (gamma_1 + gamma_2) ...
            + (1 - theta) * (gamma_1_old + gamma_2_old));
        
end
