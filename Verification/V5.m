function f = V5(DIM, i, h, h_old, phi, phi_old, k, k_old, t, PARAMS)
% Normal interior node

n = DIM.n;

% find which index we want from the re-arranged matrix
point = (DIM.r == i);
east = (DIM.r == i+1);
west = (DIM.r == i-1);
north = (DIM.r == i+n);
south = (DIM.r == i-n);

% get the dx and dz values
DELTA = DIM.DELTA(point, :);
dx = DELTA(1:2);
dz = DELTA(3:4);
ST = DIM.ST(point, :);
% get the K values
K_xx = DIM.K_xx;
K_zz = DIM.K_zz;
% get total cell volume
cell_volume = DIM.VOL(point, 5);

dt = PARAMS.dt;
theta = PARAMS.theta;
sigma = PARAMS.sigma;
source = 0;

if PARAMS.PUMPS == 1
    if DIM.XZ(point, 1) == 450 && DIM.XZ(point, 2) == 10 && phi(point) > 0.2
        source = -PARAMS.town_rate / DIM.WIDTH; % -5e-8 * 60*60*24;
    elseif DIM.XZ(point, 1) == 100 && DIM.XZ(point, 2) == 50
%         source = -PARAMS.bore_rate / DIM.WIDTH;
    end
end

% calculate k
if h(point) >= h(east)
    k_e_up   = k(point);
    k_e_down = k(east);
else
    k_e_up   = k(east);
    k_e_down = k(point);
end

if h(point) >= h(west)
    k_w_up   = k(point);
    k_w_down = k(west);
else
    k_w_up   = k(west);
    k_w_down = k(point);
end

if h(point) >= h(north) + dz(2)
    k_n_up   = k(point);
    k_n_down = k(north);
else
    k_n_up   = k(north);
    k_n_down = k(point);
end

if h(point) >= h(south) - dz(1)
    k_s_up   = k(point);
    k_s_down = k(south);
else
    k_s_up   = k(north);
    k_s_down = k(south);
end

k_e = k_e_up + sigma / 2 * (k_e_down - k_e_up);
k_n = k_n_up + sigma / 2 * (k_n_down - k_n_up);
k_w = k_w_up + sigma / 2 * (k_w_down - k_w_up);
k_s = k_s_up + sigma / 2 * (k_s_down - k_s_up);

if h_old(point) >= h_old(east)
    k_e_up_old   = k_old(point);
    k_e_down_old = k_old(east);
else
    k_e_up_old   = k_old(east);
    k_e_down_old = k_old(point);
end

if h_old(point) >= h_old(west)
    k_w_up_old   = k_old(point);
    k_w_down_old = k_old(west);
else
    k_w_up_old   = k_old(west);
    k_w_down_old = k_old(point);
end

if h_old(point) >= h_old(north) + dz(2)
    k_n_up_old   = k_old(point);
    k_n_down_old = k_old(north);
else
    k_n_up_old   = k_old(north);
    k_n_down_old = k_old(point);
end

if h_old(point) >= h_old(south) - dz(1)
    k_s_up_old   = k_old(point);
    k_s_down_old = k_old(south);
else
    k_s_up_old   = k_old(south);
    k_s_down_old = k_old(point);
end

k_e_old = k_e_up_old + sigma / 2 * (k_e_down_old - k_e_up_old);
k_n_old = k_n_up_old + sigma / 2 * (k_n_down_old - k_n_up_old);
k_w_old = k_w_up_old + sigma / 2 * (k_w_down_old - k_w_up_old);
k_s_old = k_s_up_old + sigma / 2 * (k_s_down_old - k_s_up_old);

% calculate line integrals
gamma_1 = -k_e * K_xx(ST(1)) * dz(2) / 2 * ((h(east) - h(point))/dx(2));
gamma_2 = -k_n * K_zz(ST(1)) * dx(2) / 2 * ((h(north) - h(point))/dz(2) + 1);
gamma_3 = -k_n * K_zz(ST(2)) * dx(1) / 2 * ((h(north) - h(point))/dz(2) + 1);
gamma_4 = -k_w * K_xx(ST(2)) * dz(2) / 2 * ((h(west) - h(point))/dx(1));
gamma_5 = -k_w * K_xx(ST(3)) * dz(1) / 2 * ((h(west) - h(point))/dx(1));
gamma_6 = -k_s * K_zz(ST(3)) * dx(1) / 2 * ((h(south) - h(point))/dz(1) - 1);
gamma_7 = -k_s * K_zz(ST(4)) * dx(2) / 2 * ((h(south) - h(point))/dz(1) - 1);
gamma_8 = -k_e * K_xx(ST(4)) * dz(1) / 2 * ((h(east) - h(point))/dx(2));

gamma_1_old = -k_e_old * K_xx(ST(1)) * dz(2) / 2 * ((h_old(east) - h_old(point))/dx(2));
gamma_2_old = -k_n_old * K_zz(ST(1)) * dx(2) / 2 * ((h_old(north) - h_old(point))/dz(2) + 1);
gamma_3_old = -k_n_old * K_zz(ST(2)) * dx(1) / 2 * ((h_old(north) - h_old(point))/dz(2) + 1);
gamma_4_old = -k_w_old * K_xx(ST(2)) * dz(2) / 2 * ((h_old(west) - h_old(point))/dx(1));
gamma_5_old = -k_w_old * K_xx(ST(3)) * dz(1) / 2 * ((h_old(west) - h_old(point))/dx(1));
gamma_6_old = -k_s_old * K_zz(ST(3)) * dx(1) / 2 * ((h_old(south) - h_old(point))/dz(1) - 1);
gamma_7_old = -k_s_old * K_zz(ST(4)) * dx(2) / 2 * ((h_old(south) - h_old(point))/dz(1) - 1);
gamma_8_old = -k_e_old * K_xx(ST(4)) * dz(1) / 2 * ((h_old(east) - h_old(point))/dx(2));

% evaluate f function
f = phi(point) - phi_old(point) - dt * (theta * source + (1 - theta) * source) ...
            + dt/cell_volume * (theta * ...
             (gamma_1 + gamma_2 + gamma_3 + gamma_4 + gamma_5 + gamma_6 + gamma_7 + gamma_8) ...
            + (1 - theta) * ( ...
              gamma_1_old + gamma_2_old + gamma_3_old + gamma_4_old + gamma_5_old + gamma_6_old + gamma_7_old + gamma_8_old));
        
end
