function f = V5(DIM, i, h, h_old, phi, phi_old, k, k_old, PARAMS)
% Normal interior node

n = DIM.n;

point = (DIM.r == i);
east = (DIM.r == i+1);
west = (DIM.r == i-1);
north = (DIM.r == i+n);
south = (DIM.r == i-n);

DELTA = DIM.DELTA(point, :);
dx = DELTA(1:2);
dz = DELTA(3:4);
K_xx = DIM.K_xx(point, :);
K_zz = DIM.K_zz(point, :);
cell_volume = DIM.VOL(point, 5);

dt = PARAMS.dt;
theta = PARAMS.theta;
sink = 0;

if PARAMS.PUMPS == 1
    if DIM.XZ(point, 1) == 450 && DIM.XZ(point, 2) == 10
        sink = PARAMS.town_rate * 10;
    elseif DIM.XZ(point, 1) == 100 && DIM.XZ(point, 2) == 50
        sink = PARAMS.bore_rate * 10;
    end
end

k_e = (k(point) + k(east)) / 2;
k_w = (k(point) + k(west)) / 2;
k_n = (k(point) + k(north)) / 2;
k_s = (k(point) + k(south)) / 2;

k_e_old = (k_old(point) + k_old(east)) / 2;
k_w_old = (k_old(point) + k_old(west)) / 2;
k_n_old = (k_old(point) + k_old(north)) / 2;
k_s_old = (k_old(point) + k_old(south)) / 2;

gamma_1 = k_e * K_xx(1) * dz(2) / 2 * ((h(east) - h(point))/dx(2));
gamma_2 = k_n * K_zz(1) * dx(2) / 2 * (1 + (h(north) - h(point))/dz(2));
gamma_3 = k_n * K_zz(2) * dx(1) / 2 * (1 + (h(north) - h(point))/dz(2));
gamma_4 = k_w * K_xx(2) * dz(2) / 2 * ((h(west) - h(point))/dx(1));
gamma_5 = k_w * K_xx(3) * dz(1) / 2 * ((h(west) - h(point))/dx(1));
gamma_6 = k_s * K_zz(3) * dx(1) / 2 * ((h(south) - h(point))/dz(1) - 1);
gamma_7 = k_s * K_zz(4) * dx(2) / 2 * ((h(south) - h(point))/dz(1) - 1);
gamma_8 = k_e * K_xx(4) * dz(1) / 2 * ((h(east) - h(point))/dx(2));

gamma_1_old = k_e_old * K_xx(1) * dz(2) / 2 * ((h_old(east) - h_old(point))/dx(2));
gamma_2_old = k_n_old * K_zz(1) * dx(2) / 2 * (1 + (h_old(north) - h_old(point))/dz(2));
gamma_3_old = k_n_old * K_zz(2) * dx(1) / 2 * (1 + (h_old(north) - h_old(point))/dz(2));
gamma_4_old = k_w_old * K_xx(2) * dz(2) / 2 * ((h_old(west) - h_old(point))/dx(1));
gamma_5_old = k_w_old * K_xx(3) * dz(1) / 2 * ((h_old(west) - h_old(point))/dx(1));
gamma_6_old = k_s_old * K_zz(3) * dx(1) / 2 * ((h_old(south) - h_old(point))/dz(1) - 1);
gamma_7_old = k_s_old * K_zz(4) * dx(2) / 2 * ((h_old(south) - h_old(point))/dz(1) - 1);
gamma_8_old = k_e_old * K_xx(4) * dz(1) / 2 * ((h_old(east) - h_old(point))/dx(2));

f = phi(point) - phi_old(point) - dt * (theta * sink + (1 - theta) * sink) ...
            - dt/cell_volume * (theta * ...
             (gamma_1 + gamma_2 + gamma_3 + gamma_4 + gamma_5 + gamma_6 + gamma_7 + gamma_8) ...
            + (1 - theta) * ( ...
              gamma_1_old + gamma_2_old + gamma_3_old + gamma_4_old + gamma_5_old + gamma_6_old + gamma_7_old + gamma_8_old));
        
end
