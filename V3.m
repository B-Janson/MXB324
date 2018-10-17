function f = V3(DIM, h, h_old, phi, phi_old, k, k_old, PARAMS)
% Bottom Right Corner

n = DIM.n;

point = (DIM.r == n);
west = (DIM.r == n-1);
north = (DIM.r == n+n);

DELTA = DIM.DELTA(point, :);
dx = DELTA(1);
dz = DELTA(4);
K_xx = DIM.K_xx(point, 2);
K_zz = DIM.K_zz(point, 2);
cell_volume = DIM.VOL(point, 5);

dt = PARAMS.dt;
theta = PARAMS.theta;

k_w = (k(point) + k(west)) / 2;
k_n = (k(point) + k(north)) / 2;

k_w_old = (k_old(point) + k_old(west)) / 2;
k_n_old = (k_old(point) + k_old(north)) / 2;

gamma_1 = -k_w * K_xx * dz / 2 * ((h(west) - h(point))/dx);
gamma_2 = -k_n * K_zz * dx / 2 * (1 + (h(north) - h(point))/dz);

gamma_1_old = -k_w_old * K_xx * dz / 2 * ((h_old(west) - h_old(point))/dx);
gamma_2_old = -k_n_old * K_zz * dx / 2 * (1 + (h_old(north) - h_old(point))/dz);

f = phi(point) - phi_old(point) + dt/cell_volume * (theta * ( ...
              gamma_1 + gamma_2) ...
            + (1 - theta) * ( ...
              gamma_1_old + gamma_2_old));
        
end
