function f = V7(DIM, h, h_old, phi, phi_old, k, k_old, RT, PARAMS)
% Top Left

n = DIM.n;
m=DIM.m;

% find which index we want from the re-arranged matrix
<<<<<<< HEAD
point = (DIM.r == i);
east = (DIM.r == i+1);
south = (DIM.r == i-n);
=======
point = find((DIM.r == n*(m-1)+1));
east = find((DIM.r == n*(m-1)+2));
south = find((DIM.r == n*(m-2)+1));
>>>>>>> master

% get the dx and dz values
DELTA = DIM.DELTA(point, :);
dx = DELTA(2);
dz = DELTA(3);
% get the K values
ST = DIM.ST(point, 4);
BC = DIM.BC(point, :);
K_xx = DIM.K_xx(ST);
K_zz = DIM.K_zz(ST);
% get total cell volume
cell_volume = DIM.VOL(point, 5);

dt = PARAMS.dt;
theta = PARAMS.theta;
sigma = PARAMS.sigma;

SAT = 1;

% if the water table has reached this point, don't keep pouring in water
if h_old(point) >= 0
   SAT = 0; 
end

% calculate k
if h(point) >= h(east)
    k_e_up   = k(point);
    k_e_down = k(east);
else
    k_e_up   = k(east);
    k_e_down = k(point);
end

if h(point) >= h(south) - dz
    k_s_up   = k(point);
    k_s_down = k(south);
else
    k_s_up   = k(south);
    k_s_down = k(point);
end

k_e = k_e_up + sigma / 2 * (k_e_down - k_e_up);
k_s = k_s_up + sigma / 2 * (k_s_down - k_s_up);

if h_old(point) >= h_old(east)
    k_e_up_old   = k_old(point);
    k_e_down_old = k_old(east);
else
    k_e_up_old   = k_old(east);
    k_e_down_old = k_old(point);
end

if h_old(point) >= h_old(south) - dz
    k_s_up_old   = k_old(point);
    k_s_down_old = k_old(south);
else
    k_s_up_old   = k_old(south);
    k_s_down_old = k_old(point);
end

k_e_old = k_e_up_old + sigma / 2 * (k_e_down_old - k_e_up_old);
k_s_old = k_s_up_old + sigma / 2 * (k_s_down_old - k_s_up_old);

% calculate line integrals
gamma_1 = -k_e * K_xx * dz / 2 * ((h(east) - h(point))/dx);
<<<<<<< HEAD
gamma_2 = -active_flow * r_f * dx / 2;
gamma_3 = BC(2) * dz / 2;
gamma_4 = -k_s * K_zz * dx / 2 * ((h(south) - h(point))/dz - 1);

gamma_1_old = -k_e_old * K_xx * dz / 2 * ((h_old(east) - h_old(point))/dx);
gamma_2_old = -active_flow * r_f * dx / 2;
gamma_4_old = -k_s_old * K_zz * dx / 2 * ((h_old(south) - h_old(point))/dz - 1);
=======
gamma_2 = -SAT * RT * dx / 2;
gamma_3 = -k_s * K_zz * dx / 2 * ((h(south) - h(point))/dz - 1);

gamma_1_old = -k_e_old * K_xx * dz / 2 * ((h_old(east) - h_old(point))/dx);
gamma_2_old = -SAT * RT * dx / 2;
gamma_3_old = -k_s_old * K_zz * dx / 2 * ((h_old(south) - h_old(point))/dz - 1);
>>>>>>> master

% evaluate f function
f = phi(point) - phi_old(point) + dt/cell_volume * ... 
            (theta * (gamma_1 + gamma_2 + gamma_3 + gamma_4) + ...
            (1 - theta) * (gamma_1_old + gamma_2_old + gamma_3 + gamma_4_old));
        
end
