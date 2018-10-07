function f = V1(DIM, h, h_old, phi, phi_old, k, k_old, PARAMS)
% Bottom Left Corner

% XYN = DIM.XY;
% XYN(:,DIM.r) = XYN;

n = DIM.n;
D=DIM.DELTA(1,:);
K_xx = DIM.SP(DIM.ST(1,2),2);
K_zz = DIM.SP(DIM.ST(1,1),2);
VOL = DIM.VOL(1,5);

point = (DIM.r == 1);
east = (DIM.r == 2);
north = (DIM.r == n+1);

dt = PARAMS.dt;
theta = PARAMS.theta;

gamma_1 = -(k(point) + k(east)) * K_xx * D(4) / 4 * ((h(east) - h(point))/D(2));
gamma_2 = -(k(point) + k(north)) * K_zz * D(2) / 4 * (1 + (h(north) - h(point))/D(4));

gamma_1_old = -(k_old(point) + k_old(east)) * K_xx * D(4) / 4 * ((h_old(east) - h_old(point))/D(2));
gamma_2_old = -(k_old(point) + k_old(north)) * K_zz * D(2) / 4 * (1 + (h_old(north) - h_old(point))/D(4));

f = phi(point) - phi_old(point) + dt/VOL * (theta * ( ...
              gamma_1 + gamma_2) ...
            + (1 - theta) * ( ...
              gamma_1_old + gamma_2_old));
        
end
