function f = V2(DIM, i, h, h_old, phi, phi_old, k, k_old, r_f, PARAMS)
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
% get the parameter values
ST = DIM.ST(point, :);
BC = BOUNDARY_CONDITIONS(DIM, PARAMS, DIM.XZ(point, :), r_f, h(point));
K_xx = DIM.K_xx;
K_zz = DIM.K_zz;
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

k_e = k_e_up + sigma / 2 * (k_e_down - k_e_up);
k_n = k_n_up + sigma / 2 * (k_n_down - k_n_up);
k_w = k_w_up + sigma / 2 * (k_w_down - k_w_up);

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

k_e_old = k_e_up_old + sigma / 2 * (k_e_down_old - k_e_up_old);
k_n_old = k_n_up_old + sigma / 2 * (k_n_down_old - k_n_up_old);
k_w_old = k_w_up_old + sigma / 2 * (k_w_down_old - k_w_up_old);

% calculate line integrals
gamma_1 = -k_e * K_xx(ST(1)) * dz(2) / 2 * ((h(east) - h(point))/dx(2));
gamma_2 = -k_n * K_zz(ST(1)) * dx(2) / 2 * ((h(north) - h(point))/dz(2) + 1);
gamma_3 = -k_n * K_zz(ST(2)) * dx(1) / 2 * ((h(north) - h(point))/dz(2) + 1);
gamma_4 = -k_w * K_xx(ST(2)) * dz(2) / 2 * ((h(west) - h(point))/dx(1));
gamma_5 = BC(1) * dx(1) / 2;
gamma_6 = BC(2) * dx(2) / 2;

gamma_1_old = -k_e_old * K_xx(ST(1)) * dz(2) / 2 * ((h_old(east) - h_old(point))/dx(2));
gamma_2_old = -k_n_old * K_zz(ST(1)) * dx(2) / 2 * ((h_old(north) - h_old(point))/dz(2) + 1);
gamma_3_old = -k_n_old * K_zz(ST(2)) * dx(1) / 2 * ((h_old(north) - h_old(point))/dz(2) + 1);
gamma_4_old = -k_w_old * K_xx(ST(2)) * dz(2) / 2 * ((h_old(west) - h_old(point))/dx(1));

% evaluate f function
f = phi(point) - phi_old(point) + dt/cell_volume * (theta * ( ...
              gamma_1 + gamma_2 + gamma_3 + gamma_4 + gamma_5 + gamma_6) ...
            + (1 - theta) * ( ...
              gamma_1_old + gamma_2_old + gamma_3_old + gamma_4_old + gamma_5 + gamma_6));
        
end
