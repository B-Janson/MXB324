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

f = phi(point) - phi_old(point) - dt/cell_volume * (theta * ( ...
              (k(point) + k(west)) * K_xx * dz/(4 * dx) * (h(west) - h(point)) ...
            + (k(point) + k(south)) * K_zz * dx / 4 * ((h(south) - h(point)) / dz - 1)  ...
            + (r_f + r_f * cos(2*pi*t/365)) * dx/2) ...
            + (1 - theta) * ( ...
              (k_old(point) + k_old(west)) * K_xx * dz/(4 * dx) * (h_old(west) - h_old(point)) ...
            + (k_old(point) + k_old(south)) * K_zz * dx / 4 * ((h_old(south) - h_old(point)) / dz - 1)  ...
            + (r_f + r_f * cos(2*pi*(t-dt)/365)) * dx/2));
        
end
