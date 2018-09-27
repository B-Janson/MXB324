function f = V1(DIM, h, h_old, phi, phi_old, k, k_old, dt, theta) %, sigma)
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

f = phi(point) - phi_old(point) - dt/cell_volume * (theta * ( ...
              (k(point) + k(east)) * K_xx * dz/(4 * dx) * (h(east) - h(point)) ...
            + (k(point) + k(north)) * K_zz * dx / 4 * (1 + (h(north) - h(point)) / dz))  ...
            + (1 - theta) * ( ...
              (k_old(point) + k_old(east)) * K_xx * dz/(4 * dx) * (h_old(east) - h_old(point)) ...
            + (k_old(point) + k_old(north)) * K_zz * dx / 4 * (1 + (h_old(north) - h_old(point)) / dz)));
        
end
