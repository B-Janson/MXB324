function f = V29(DIM, i, h, h_old, phi, phi_old, k, k_old, dt, t, theta, r_f) %, sigma)
% Sandstone/Rain

n = DIM.n;
dx = DIM.dx(1,1);
dz = DIM.dz(1,1);
K_xx = DIM.K_xx(1,1);
K_zz = DIM.K_zz(1,1);
cell_volume = DIM.cell_volume(1,1);

point = (DIM.r == i);
east = (DIM.r == i+1);
south = (DIM.r == i-n);

f = phi(point) - phi_old(point) - dt/cell_volume * (theta * ( ...
              (k(point) + k(east)) * K_xx * dz/(4 * dx) * (h(east) - h(point)) ...
            + (k(point) + k(south)) * K_zz * dx / 4 * ((h(south) - h(point)) / dz - 1)  ...
            + (r_f + r_f * cos(2*pi*t/365)) * dx/2) ...
            + (1 - theta) * ( ...
              (k_old(point) + k_old(east)) * K_xx * dz/(4 * dx) * (h_old(east) - h_old(point)) ...
            + (k_old(point) + k_old(south)) * K_zz * dx / 4 * ((h_old(south) - h_old(point)) / dz - 1)  ...
            + (r_f + r_f * cos(2*pi*(t-dt)/365)) * dx/2));
        
end
