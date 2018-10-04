function f = V1(DIM, h, h_old, phi, phi_old, k, k_old, PARAMS)
% Bottom Left Corner

% XYN = DIM.XY;
% XYN(:,DIM.r) = XYN;

n = DIM.n;
dx = DIM.dx(1,1);
dz = DIM.dz(1,1);
K_xx = DIM.K_xx(1,1);
K_zz = DIM.K_zz(1,1);
cell_volume = DIM.cell_volume(1,1);

point = (DIM.r == 1);
east = (DIM.r == 2);
north = (DIM.r == n+1);

dt = PARAMS.dt;
theta = PARAMS.theta;

gamma_1 = -(k(point) + k(east)) * K_xx * dz / 4 * ((h(east) - h(point))/dx);
gamma_2 = -(k(point) + k(north)) * K_zz * dx / 4 * (1 + (h(north) - h(point))/dz);

gamma_1_old = -(k_old(point) + k_old(east)) * K_xx * dz / 4 * ((h_old(east) - h_old(point))/dx);
gamma_2_old = -(k_old(point) + k_old(north)) * K_zz * dx / 4 * (1 + (h_old(north) - h_old(point))/dz);

f = phi(point) - phi_old(point) + dt/cell_volume * (theta * ( ...
              gamma_1 + gamma_2) ...
            + (1 - theta) * ( ...
              gamma_1_old + gamma_2_old));
        
end
