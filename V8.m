function f = V8(DIM, i, h, h_old, phi, phi_old, k, k_old, PARAMS)
% Right boundary sandstone

% XYN = DIM.XY;
% XYN(:,DIM.r) = XYN;

n = DIM.n;
D=DIM.DELTA(i,:);
K_xx = DIM.SP(DIM.ST(i,2),2);
K_zz = DIM.SP(DIM.ST(i,1),2);
VOL = DIM.VOL(i,5);

point = (DIM.r == i);
west = (DIM.r == i-1);
north = (DIM.r == i+n);
south = (DIM.r == i-n);

dt = PARAMS.dt;
theta = PARAMS.theta;

f = phi(point) - phi_old(point) - dt/VOL * (theta * ( ...
              (k(point) + k(west)) * K_xx * D(4)/(4 * D(1)) * (h(west) - h(point)) ...
            + (k(point) + k(north)) * K_zz * D(1) / 4 * ((h(north) - h(point)) / D(4) + 1)  ...
            + (k(point) + k(south)) * K_zz * D(1) / 4 * ((h(south) - h(point)) / D(3) - 1)  ...
            + (k(point) + k(west)) * K_xx * D(4)/(4 * D(1)) * (h(west) - h(point))) ...
            + (1 - theta) * ( ...
              (k_old(point) + k_old(west)) * K_xx * D(4)/(4 * D(1)) * (h_old(west) - h_old(point)) ...
            + (k_old(point) + k_old(north)) * K_zz * D(1) / 4 * ((h_old(north) - h_old(point)) / D(4) + 1)  ...
            + (k_old(point) + k_old(south)) * K_zz * D(1) / 4 * ((h_old(south) - h_old(point)) / D(3) - 1)  ...
            + (k_old(point) + k_old(west)) * K_xx * D(4)/(4 * D(1)) * (h_old(west) - h_old(point))));
        
end
