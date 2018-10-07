function f = V3(DIM, h, h_old, phi, phi_old, k, k_old, PARAMS)
% Bottom Right Corner

% XYN = DIM.XY;
% XYN(:,DIM.r) = XYN;

n = DIM.n;
dx = DIM.dx(1,1);
dz = DIM.dz(1,1);
K_xx = DIM.K_xx(1,1);
K_zz = DIM.K_zz(1,1);
VOL = DIM.VOL(1);
point = (DIM.r == n);
west = (DIM.r == n-1);
north = (DIM.r == n+n);

dt = PARAMS.dt;
theta = PARAMS.theta;

gamma_1 = -(k(point) + k(north)) * K_zz * dx * (1 + (h(north) - h(point))/dz) / 4;
gamma_2 = -(k(point) + k(west)) * K_xx * dz * ((h(west) - h(point))/dx) / 4;

gamma_1_old = -(k_old(point) + k_old(north)) * K_zz * dx * (1 + (h_old(north) - h_old(point))/dz) / 4;
gamma_2_old = -(k_old(point) + k_old(west)) * K_xx * dz * ((h_old(west) - h_old(point))/dx) / 4;

f = phi(point) - phi_old(point) + dt/VOL * (theta * ( ...
              gamma_1 + gamma_2) ...
            + (1 - theta) * ( ...
              gamma_1_old + gamma_2_old));
        
end
