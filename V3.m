function f = V3(DIM, h, h_old, phi, phi_old, k, k_old, PARAMS)
% Bottom Right Corner

% XYN = DIM.XY;
% XYN(:,DIM.r) = XYN;

n = DIM.n;
D=DIM.DELTA(n,:);
K_xx = DIM.SP(DIM.ST(n,2),2);
K_zz = DIM.SP(DIM.ST(n,1),2);
VOL = DIM.VOL(n,5);

point = (DIM.r == n);
west = (DIM.r == n-1);
north = (DIM.r == n+n);

dt = PARAMS.dt;
theta = PARAMS.theta;

gamma_1 = -(k(point) + k(north)) * K_zz * D(1) * (1 + (h(north) - h(point))/D(4)) / 4;
gamma_2 = -(k(point) + k(west)) * K_xx * D(4) * ((h(west) - h(point))/D(1)) / 4;

gamma_1_old = -(k_old(point) + k_old(north)) * K_zz * D(1) * (1 + (h_old(north) - h_old(point))/D(4)) / 4;
gamma_2_old = -(k_old(point) + k_old(west)) * K_xx * D(4) * ((h_old(west) - h_old(point))/D(1)) / 4;

f = phi(point) - phi_old(point) + dt/VOL * (theta * ( ...
              gamma_1 + gamma_2) ...
            + (1 - theta) * ( ...
              gamma_1_old + gamma_2_old));
        
end
