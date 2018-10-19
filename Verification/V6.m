function f = V6(DIM, i, h, h_old, phi, phi_old, k, k_old, PARAMS)
% Right boundary

n = DIM.n;

% find which index we want from the re-arranged matrix
point = (DIM.r == i);
west = (DIM.r == i-1);
north = (DIM.r == i+n);
south = (DIM.r == i-n);

% get the dx and dz values
DELTA = DIM.DELTA(point, :);
dx = DELTA(1);
dz = DELTA(3:4);
% get the K values
ST = DIM.ST(point, :);
K_xx = DIM.K_xx;
K_zz = DIM.K_zz;
% get total cell volume
cell_volume = DIM.VOL(point, 5);

dt = PARAMS.dt;
theta = PARAMS.theta;
sigma = PARAMS.sigma;

% calculate k
k_n_up   = max(k(point), k(north));
k_n_down = min(k(point), k(north));
k_w_up   = max(k(point), k(west));
k_w_down = min(k(point), k(west));
k_s_up   = max(k(point), k(south));
k_s_down = min(k(point), k(south));

k_n = k_n_up + sigma / 2 * (k_n_down - k_n_up);
k_w = k_w_up + sigma / 2 * (k_w_down - k_w_up);
k_s = k_s_up + sigma / 2 * (k_s_down - k_s_up);

k_n_up_old   = max(k_old(point), k_old(north));
k_n_down_old = min(k_old(point), k_old(north));
k_w_up_old   = max(k_old(point), k_old(west));
k_w_down_old = min(k_old(point), k_old(west));
k_s_up_old   = max(k_old(point), k_old(south));
k_s_down_old = min(k_old(point), k_old(south));

k_n_old = k_n_up_old + sigma / 2 * (k_n_down_old - k_n_up_old);
k_w_old = k_w_up_old + sigma / 2 * (k_w_down_old - k_w_up_old);
k_s_old = k_s_up_old + sigma / 2 * (k_s_down_old - k_s_up_old);

% calculate line integrals
gamma_1 = -k_n * K_zz(ST(2)) * dx / 2 * ((h(north) - h(point))/dz(2) + 1);
gamma_2 = -k_w * K_xx(ST(2)) * dz(2) / 2 * ((h(west) - h(point))/dx);
gamma_3 = -k_w * K_xx(ST(3)) * dz(1) / 2 * ((h(west) - h(point))/dx);
gamma_4 = -k_s * K_zz(ST(3)) * dx / 2 * ((h(south) - h(point))/dz(1) - 1);

gamma_1_old = -k_n_old * K_zz(ST(2)) * dx / 2 * ((h_old(north) - h_old(point))/dz(2) + 1);
gamma_2_old = -k_w_old * K_xx(ST(2)) * dz(2) / 2 * ((h_old(west) - h_old(point))/dx);
gamma_3_old = -k_w_old * K_xx(ST(3)) * dz(1) / 2 * ((h_old(west) - h_old(point))/dx);
gamma_4_old = -k_s_old * K_zz(ST(3)) * dx / 2 * ((h_old(south) - h_old(point))/dz(1) - 1);

% evaluate f function
f = phi(point) - phi_old(point) + dt/cell_volume * (theta * ...
             (gamma_1 + gamma_2 + gamma_3 + gamma_4) ...
            + (1 - theta) * ...
              (gamma_1_old + gamma_2_old + gamma_3_old + gamma_4_old));
        
end
