function f = V2(DIM, i, h, h_old, phi, phi_old, k, k_old, PARAMS)
% V2  returns the f function evaulated at the bottom horizontal boundary of the
%     grid

n = DIM.n;

% find which index we want from the re-arranged matrix
point = (DIM.r == i);
east = (DIM.r == i+1);
west = (DIM.r == i-1);
north = (DIM.r == i+n);

% get the dx and dz values
DELTA = DIM.DELTA(point, :);
dx = DELTA(1:2);
dz = DELTA(3:4);
% get the K values
K_xx = DIM.K_xx(point, :);
K_zz = DIM.K_zz(point, :);
% get total cell volume
cell_volume = DIM.VOL(point, 5);

dt = PARAMS.dt;
theta = PARAMS.theta;

% calculate k
k_e = (k(point) + k(east)) / 2;
k_w = (k(point) + k(west)) / 2;
k_n = (k(point) + k(north)) / 2;

k_e_old = (k_old(point) + k_old(east)) / 2;
k_w_old = (k_old(point) + k_old(west)) / 2;
k_n_old = (k_old(point) + k_old(north)) / 2;

% calculate line integrals
gamma_1 = -k_e * K_xx(1) * dz(2) / 2 * ((h(east) - h(point))/dx(2));
gamma_2 = -k_n * K_zz(1) * dx(2) / 2 * (1 + (h(north) - h(point))/dz(2));
gamma_3 = -k_n * K_zz(2) * dx(1) / 2 * (1 + (h(north) - h(point))/dz(2));
gamma_4 = -k_w * K_xx(2) * dz(2) / 2 * ((h(west) - h(point))/dx(1));

gamma_1_old = -k_e_old * K_xx(1) * dz(2) / 2 * ((h_old(east) - h_old(point))/dx(2));
gamma_2_old = -k_n_old * K_zz(1) * dx(2) / 2 * (1 + (h_old(north) - h_old(point))/dz(2));
gamma_3_old = -k_n_old * K_zz(2) * dx(1) / 2 * (1 + (h_old(north) - h_old(point))/dz(2));
gamma_4_old = -k_w_old * K_xx(2) * dz(2) / 2 * ((h_old(west) - h_old(point))/dx(1));

% evaluate f function
f = phi(point) - phi_old(point) + dt/cell_volume * (theta * ( ...
              gamma_1 + gamma_2 + gamma_3 + gamma_4) ...
            + (1 - theta) * ( ...
              gamma_1_old + gamma_2_old + gamma_3_old + gamma_4_old));
        
end
