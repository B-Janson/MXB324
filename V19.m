function f = V19(DIM, i, h, h_old, phi, phi_old, k, k_old, dt, theta) %, sigma)
% Normal sandstone interior node

n = DIM.n;
dx = DIM.dx(1,1);
dz = DIM.dz(1,1);
K_xx = DIM.K_xx(1,1);
K_zz = DIM.K_zz(1,1);
cell_volume = sum(DIM.cell_volume(i,:));

point = (DIM.r == i);
east = (DIM.r == i+1);
west = (DIM.r == i-1);
north = (DIM.r == i+n);
south = (DIM.r == i-n);

f = phi(point) - phi_old(point) - dt/cell_volume * (theta * ( ...
              (k(point) + k(east)) * K_xx * dz/(4 * dx) * (h(east) - h(point)) ... % gamma 1
            + (k(point) + k(north)) * K_zz * dx / 4 * ((h(north) - h(point)) / dz + 1)  ... % gamma 2
            + (k(point) + k(north)) * K_zz * dx / 4 * ((h(north) - h(point)) / dz + 1)  ... % gamma 3
            + (k(point) + k(west)) * K_xx * dz/(4 * dx) * (h(west) - h(point)) ... % gamma 4
            + (k(point) + k(west)) * K_xx * dz/(4 * dx) * (h(west) - h(point)) ... % gamma 5
            + (k(point) + k(south)) * K_zz * dx / 4 * ((h(south) - h(point)) / dz - 1)  ... % gamma 6
            + (k(point) + k(south)) * K_zz * dx / 4 * ((h(south) - h(point)) / dz - 1)  ... % gamma 7
            + (k(point) + k(east)) * K_xx * dz/(4 * dx) * (h(east) - h(point))) ... % gamma 8
            + (1 - theta) * ( ...
              (k_old(point) + k_old(east)) * K_xx * dz/(4 * dx) * (h_old(east) - h_old(point)) ... % gamma 1
            + (k_old(point) + k_old(north)) * K_zz * dx / 4 * ((h_old(north) - h_old(point)) / dz + 1)  ... % gamma 2
            + (k_old(point) + k_old(north)) * K_zz * dx / 4 * ((h_old(north) - h_old(point)) / dz + 1)  ... % gamma 3
            + (k_old(point) + k_old(west)) * K_xx * dz/(4 * dx) * (h_old(west) - h_old(point)) ... % gamma 4
            + (k_old(point) + k_old(west)) * K_xx * dz/(4 * dx) * (h_old(west) - h_old(point)) ... % gamma 5
            + (k_old(point) + k_old(south)) * K_zz * dx / 4 * ((h_old(south) - h_old(point)) / dz - 1)  ... % gamma 6
            + (k_old(point) + k_old(south)) * K_zz * dx / 4 * ((h_old(south) - h_old(point)) / dz - 1)  ... % gamma 7
            + (k_old(point) + k_old(east)) * K_xx * dz/(4 * dx) * (h_old(east) - h_old(point)))); % gamma 8
        
end
