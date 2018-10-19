function f = V3(DIM, h, h_old, phi, phi_old, k, k_old, PARAMS)
% V3  returns the f function evaulated at the bottom right corner of the
%     grid

n = DIM.n;

% find which index we want from the re-arranged matrix
point = (DIM.r == n);
west = (DIM.r == n-1);
north = (DIM.r == n+n);

% get the dx and dz values
DELTA = DIM.DELTA(point, :);
dx = DELTA(1);
dz = DELTA(4);
% get the K values for the second quadrant only
ST = DIM.ST(point, 2);
K_xx = DIM.K_xx(ST);
K_zz = DIM.K_zz(ST);
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

k_n = k_n_up + sigma / 2 * (k_n_down - k_n_up);
k_w = k_w_up + sigma / 2 * (k_w_down - k_w_up);

k_n_up_old   = max(k_old(point), k_old(north));
k_n_down_old = min(k_old(point), k_old(north));
k_w_up_old   = max(k_old(point), k_old(west));
k_w_down_old = min(k_old(point), k_old(west));

k_n_old = k_n_up_old + sigma / 2 * (k_n_down_old - k_n_up_old);
k_w_old = k_w_up_old + sigma / 2 * (k_w_down_old - k_w_up_old);

gamma_1 = -k_w * K_xx * dz / 2 * ((h(west) - h(point))/dx);
gamma_2 = -k_n * K_zz * dx / 2 * ((h(north) - h(point))/dz + 1);

% calculate line integrals
gamma_1_old = -k_w_old * K_xx * dz / 2 * ((h_old(west) - h_old(point))/dx);
gamma_2_old = -k_n_old * K_zz * dx / 2 * ((h_old(north) - h_old(point))/dz + 1);

% evaluate f function
f = phi(point) - phi_old(point) + dt/cell_volume * (theta * ( ...
              gamma_1 + gamma_2) ...
            + (1 - theta) * ( ...
              gamma_1_old + gamma_2_old));
        
end
