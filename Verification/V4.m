function f = V4(DIM, i, h, h_old, phi, phi_old, k, k_old, PARAMS)
% V4  returns the f function evaulated at the left boundary of the
%     grid

n = DIM.n;

% find which index we want from the re-arranged matrix
point = (DIM.r == i);
east = (DIM.r == i+1);
north = (DIM.r == i+n);
south = (DIM.r == i-n);

% get the dx and dz values
DELTA = DIM.DELTA(point, :);
dx = DELTA(2);
dz = DELTA(3:4);

ST = DIM.ST(point, :);
K_xx = DIM.K_xx;
K_zz = DIM.K_zz;
% get total cell volume
cell_volume = DIM.VOL(point, 5);

dt = PARAMS.dt;
theta = PARAMS.theta;

% calculate k
k_e = (k(point) + k(east)) / 2;
k_n = (k(point) + k(north)) / 2;
k_s = (k(point) + k(south)) / 2;

k_e_old = (k_old(point) + k_old(east)) / 2;
k_n_old = (k_old(point) + k_old(north)) / 2;
k_s_old = (k_old(point) + k_old(south)) / 2;

% calculate line integrals
gamma_1 = -k_e * K_xx(ST(1)) * dz(2) / 2 * ((h(east) - h(point))/dx);
gamma_2 = -k_n * K_zz(ST(1)) * dx / 2 * ((h(north) - h(point))/dz(2) + 1);
gamma_3 = -k_s * K_zz(ST(4)) * dx / 2 * ((h(south) - h(point))/dz(1) - 1);
gamma_4 = -k_e * K_xx(ST(4)) * dz(1) / 2 * ((h(east) - h(point))/dx);

gamma_1_old = -k_e_old * K_xx(ST(1)) * dz(2) / 2 * ((h_old(east) - h_old(point))/dx);
gamma_2_old = -k_n_old * K_zz(ST(1)) * dx / 2 * ((h_old(north) - h_old(point))/dz(2) + 1);
gamma_3_old = -k_s_old * K_zz(ST(4)) * dx / 2 * ((h_old(south) - h_old(point))/dz(1) - 1);
gamma_4_old = -k_e_old * K_xx(ST(4)) * dz(1) / 2 * ((h_old(east) - h_old(point))/dx);

% evaluate f function
f = phi(point) - phi_old(point) + dt/cell_volume * ( ...
            theta * (gamma_1 + gamma_2 + gamma_3 + gamma_4) ...
            + (1 - theta) * (gamma_1_old + gamma_2_old + gamma_3_old + gamma_4_old));
        
end
