function f = V1(DIM, h, h_old, phi, phi_old, k, k_old, r_f, PARAMS)
% V1  returns the f function evaulated at the bottom left corner of the
%     grid

n = DIM.n;

% find which index we want from the re-arranged matrix
point = (DIM.r == 1);
east = (DIM.r == 2);
north = (DIM.r == n+1);

% get the dx and dz values
DELTA = DIM.DELTA(point, :);
dx = DELTA(2);
dz = DELTA(4);
% get the K values for the first quadrant only
ST = DIM.ST(point, 1);
BC = BOUNDARY_CONDITIONS(DIM, PARAMS, DIM.XZ(point, :), r_f, h(point));
K_xx = DIM.K_xx(ST);
K_zz = DIM.K_zz(ST);
% get total cell volume
cell_volume = DIM.VOL(point, 5);

dt = PARAMS.dt;
theta = PARAMS.theta;
sigma = PARAMS.sigma;

% calculate k
if h(point) >= h(east)
    k_e_up   = k(point);
    k_e_down = k(east);
else
    k_e_up   = k(east);
    k_e_down = k(point);
end

if h(point) >= h(north) + dz
    k_n_up   = k(point);
    k_n_down = k(north);
else
    k_n_up   = k(north);
    k_n_down = k(point);
end

k_e = k_e_up + sigma / 2 * (k_e_down - k_e_up);
k_n = k_n_up + sigma / 2 * (k_n_down - k_n_up);

if h_old(point) >= h_old(east)
    k_e_up_old   = k_old(point);
    k_e_down_old = k_old(east);
else
    k_e_up_old   = k_old(east);
    k_e_down_old = k_old(point);
end

if h_old(point) >= h_old(north) + dz
    k_n_up_old   = k_old(point);
    k_n_down_old = k_old(north);
else
    k_n_up_old   = k_old(north);
    k_n_down_old = k_old(point);
end

k_e_old = k_e_up_old + sigma / 2 * (k_e_down_old - k_e_up_old);
k_n_old = k_n_up_old + sigma / 2 * (k_n_down_old - k_n_up_old);

% calculate line integrals
gamma_1 = -k_e * K_xx * dz / 2 * ((h(east) - h(point))/dx);
gamma_2 = -k_n * K_zz * dx / 2 * ((h(north) - h(point))/dz + 1);
gamma_3 = BC(1) * dz / 2;
gamma_4 = BC(2) * dx / 2;

gamma_1_old = -k_e_old * K_xx * dz / 2 * ((h_old(east) - h_old(point))/dx);
gamma_2_old = -k_n_old * K_zz * dx / 2 * ((h_old(north) - h_old(point))/dz + 1);

% evaluate f function
f = phi(point) - phi_old(point) + dt/cell_volume * ( ...
              theta * (gamma_1 + gamma_2 + gamma_3 + gamma_4) ...
            + (1 - theta) * (gamma_1_old + gamma_2_old + gamma_3 + gamma_4));
        
end
