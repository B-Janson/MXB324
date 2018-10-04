function f = V2(DIM, i, h, h_old, phi, phi_old, k, k_old, PARAMS)
% Bottom Boundary

n = DIM.n;
dx = DIM.dx(1,1);
dz = DIM.dz(1,1);
K_xx = DIM.K_xx(1,1);
K_zz = DIM.K_zz(1,1);
cell_volume = DIM.cell_volume(1,1) + DIM.cell_volume(1,1);

point = (DIM.r == i);
east = (DIM.r == i+1);
west = (DIM.r == i-1);
north = (DIM.r == i+n);

dt = PARAMS.dt;
theta = PARAMS.theta;

f = phi(point) - phi_old(point) - dt/cell_volume * (theta * ( ...
              (k(point) + k(east)) * K_xx * dz/(4 * dx) * (h(east) - h(point)) ...
            + (k(point) + k(north)) * K_zz * dx / 4 * (1 + (h(north) - h(point)) / dz)  ...
            + (k(point) + k(north)) * K_zz * dx / 4 * (1 + (h(north) - h(point)) / dz)  ...
            + (k(point) + k(west)) * K_xx * dz/(4 * dx) * (h(west) - h(point))) ...
            + (1 - theta) * ( ...
              (k_old(point) + k_old(east)) * K_xx * dz/(4 * dx) * (h_old(east) - h_old(point)) ...
            + (k_old(point) + k_old(north)) * K_zz * dx / 4 * (1 + (h_old(north) - h_old(point)) / dz)  ...
            + (k_old(point) + k_old(north)) * K_zz * dx / 4 * (1 + (h_old(north) - h_old(point)) / dz)  ...
            + (k_old(point) + k_old(west)) * K_xx * dz/(4 * dx) * (h_old(west) - h_old(point))));
        
end
