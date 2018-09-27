function f = V8(DIM, i, h, h_old, phi, phi_old, k, k_old, dt, theta) %, sigma)
% Right boundary sandstone

% XYN = DIM.XY;
% XYN(:,DIM.r) = XYN;

n = DIM.n;
dx = DIM.dx(1,1);
dz = DIM.dz(1,1);
K_xx = DIM.K_xx(1,1);
K_zz = DIM.K_zz(1,1);
cell_volume = DIM.cell_volume(1,1) + DIM.cell_volume(1,1);

point = (DIM.r == i);
west = (DIM.r == i-1);
north = (DIM.r == i+n);
south = (DIM.r == i-n);

f = phi(point) - phi_old(point) - dt/cell_volume * (theta * ( ...
              (k(point) + k(west)) * K_xx * dz/(4 * dx) * (h(west) - h(point)) ...
            + (k(point) + k(north)) * K_zz * dx / 4 * ((h(north) - h(point)) / dz + 1)  ...
            + (k(point) + k(south)) * K_zz * dx / 4 * ((h(south) - h(point)) / dz - 1)  ...
            + (k(point) + k(west)) * K_xx * dz/(4 * dx) * (h(west) - h(point))) ...
            + (1 - theta) * ( ...
              (k_old(point) + k_old(west)) * K_xx * dz/(4 * dx) * (h_old(west) - h_old(point)) ...
            + (k_old(point) + k_old(north)) * K_zz * dx / 4 * ((h_old(north) - h_old(point)) / dz + 1)  ...
            + (k_old(point) + k_old(south)) * K_zz * dx / 4 * ((h_old(south) - h_old(point)) / dz - 1)  ...
            + (k_old(point) + k_old(west)) * K_xx * dz/(4 * dx) * (h_old(west) - h_old(point))));
        
end
