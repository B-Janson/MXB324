function f = V4(DIM, i, h, h_old, phi, phi_old, k, k_old, PARAMS)
% V4  returns the f function evaulated at the left boundary of the
%     grid

n = DIM.n;m=DIM.m;


% find which index we want from the re-arranged matrix
point = find((DIM.r == i));
east = find(DIM.r == i+1);
north = find(DIM.r == i+n);
south = find((DIM.r == i-n));

% get the dx and dz values
DELTA = DIM.DELTA(point, :);
dx = DELTA(2);
dz = DELTA(3:4);
% get the K values
ST = DIM.ST(point, :);
BC = DIM.BC(point, :);
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
    k_s_up   = k(south);
    k_s_down = k(point);
end

k_e = k_e_up + sigma / 2 * (k_e_down - k_e_up);
k_n = k_n_up + sigma / 2 * (k_n_down - k_n_up);
k_s = k_s_up + sigma / 2 * (k_s_down - k_s_up);

if h_old(point) >= h_old(east)
    k_e_up_old   = k_old(point);
    k_e_down_old = k_old(east);
else
    k_e_up_old   = k_old(east);
    k_e_down_old = k_old(point);
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
k_s_old = k_s_up_old + sigma / 2 * (k_s_down_old - k_s_up_old);

% calculate line integrals
gamma_1 = -k_e * K_xx(ST(1)) * dz(2) / 2 * ((h(east) - h(point))/dx);
gamma_2 = -k_n * K_zz(ST(1)) * dx / 2 * ((h(north) - h(point))/dz(2) + 1);
gamma_3 = BC(2) * dz(2) / 2;
gamma_4 = BC(1) * dz(1) / 2;
gamma_5 = -k_s * K_zz(ST(4)) * dx / 2 * ((h(south) - h(point))/dz(1) - 1);
gamma_6 = -k_e * K_xx(ST(4)) * dz(1) / 2 * ((h(east) - h(point))/dx);

gamma_1_old = -k_e_old * K_xx(ST(1)) * dz(2) / 2 * ((h_old(east) - h_old(point))/dx);
gamma_2_old = -k_n_old * K_zz(ST(1)) * dx / 2 * ((h_old(north) - h_old(point))/dz(2) + 1);
gamma_5_old = -k_s_old * K_zz(ST(4)) * dx / 2 * ((h_old(south) - h_old(point))/dz(1) - 1);
gamma_6_old = -k_e_old * K_xx(ST(4)) * dz(1) / 2 * ((h_old(east) - h_old(point))/dx);

% evaluate f function
f = phi(point) - phi_old(point) + dt/cell_volume * ( ...
            theta * (gamma_1 + gamma_2 + gamma_3 + gamma_4 + gamma_5 + gamma_6) ...
            + (1 - theta) * (gamma_1_old + gamma_2_old + gamma_3 + gamma_4 + gamma_5_old + gamma_6_old));
        
end
