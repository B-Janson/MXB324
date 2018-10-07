function f = V15(DIM, i, h, h_old, phi, phi_old, k, k_old, PARAMS)
% Left boundary alluvium

% XYN = DIM.XY;
% XYN(:,DIM.r) = XYN;

n = DIM.n;
dx = DIM.dx(1,1);
dz = DIM.dz(1,1);
K_xx = DIM.K_xx(1,1);
K_zz = DIM.K_zz(1,1);
VOL = DIM.VOL(1,1) + DIM.VOL(1,1);

point = (DIM.r == i);
east = (DIM.r == i+1);
north = (DIM.r == i+n);
south = (DIM.r == i-n);

dt = PARAMS.dt;
theta = PARAMS.theta;

f = phi(point) - phi_old(point) - dt/VOL * (theta * ( ...
              (k(point) + k(east)) * K_xx * dz/(4 * dx) * (h(east) - h(point)) ...
            + (k(point) + k(north)) * K_zz * dx / 4 * ((h(north) - h(point)) / dz + 1)  ...
            + (k(point) + k(south)) * K_zz * dx / 4 * ((h(south) - h(point)) / dz - 1)  ...
            + (k(point) + k(east)) * K_xx * dz/(4 * dx) * (h(east) - h(point))) ...
            + (1 - theta) * ( ...
              (k_old(point) + k_old(east)) * K_xx * dz/(4 * dx) * (h_old(east) - h_old(point)) ...
            + (k_old(point) + k_old(north)) * K_zz * dx / 4 * ((h_old(north) - h_old(point)) / dz + 1)  ...
            + (k_old(point) + k_old(south)) * K_zz * dx / 4 * ((h_old(south) - h_old(point)) / dz - 1)  ...
            + (k_old(point) + k_old(east)) * K_xx * dz/(4 * dx) * (h_old(east) - h_old(point))));
        
end
