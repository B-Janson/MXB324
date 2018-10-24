function f = V8(DIM, i, h, h_old, phi, phi_old, k, k_old, t, PARAMS)
% Top row

n = DIM.n;

% find which index we want from the re-arranged matrix
point = (DIM.r == i);
east = (DIM.r == i+1);
west = (DIM.r == i-1);
south = (DIM.r == i-n);

% get the dx and dz values
DELTA = DIM.DELTA(point, :);
dx = DELTA(1:2);
dz = DELTA(3);

% get the K values for the first quadrant only
ST = DIM.ST(point, :);
BC = BOUNDARY_CONDITIONS(DIM, PARAMS, DIM.XZ(point, :), t, h(point));
K_xx = DIM.K_xx;
K_zz = DIM.K_zz;
% get total cell volume
cell_volume = DIM.VOL(point, 5);

dt = PARAMS.dt;
theta = PARAMS.theta;
sigma = PARAMS.sigma;
active_flow = 1;

% if the water table has reached this point, don't keep pouring in water
if h_old(point) >= 0
   active_flow = 0; 
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

if h(point) >= h(south) - dz(1)
    k_s_up   = k(point);
    k_s_down = k(south);
else
    k_s_up   = k(south);
    k_s_down = k(point);
end

k_e = k_e_up + sigma / 2 * (k_e_down - k_e_up);
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

if h_old(point) >= h_old(south) - dz(1)
    k_s_up_old   = k_old(point);
    k_s_down_old = k_old(south);
else
    k_s_up_old   = k_old(south);
    k_s_down_old = k_old(point);
end

k_e_old = k_e_up_old + sigma / 2 * (k_e_down_old - k_e_up_old);
k_w_old = k_w_up_old + sigma / 2 * (k_w_down_old - k_w_up_old);
k_s_old = k_s_up_old + sigma / 2 * (k_s_down_old - k_s_up_old);

% calculate line integrals
gamma_1 = -k_e * K_xx(ST(4)) * dz / 2 * ((h(east) - h(point))/dx(2));
gamma_2 = active_flow * BC(2) * dx(2) / 2;
gamma_3 = active_flow * BC(1) * dx(1) / 2;
gamma_4 = -k_w * K_xx(ST(3)) * dz / 2 * ((h(west) - h(point))/dx(1));
gamma_5 = -k_s * K_zz(ST(3)) * dx(1) / 2 * ((h(south) - h(point))/dz - 1);
gamma_6 = -k_s * K_zz(ST(4)) * dx(2) / 2 * ((h(south) - h(point))/dz - 1);

gamma_1_old = -k_e_old * K_xx(ST(4)) * dz / 2 * ((h_old(east) - h_old(point))/dx(2));
gamma_2_old = active_flow * BC(2) * dx(2) / 2;
gamma_3_old = active_flow * BC(1) * dx(1) / 2;
gamma_4_old = -k_w_old * K_xx(ST(3)) * dz / 2 * ((h_old(west) - h_old(point))/dx(1));
gamma_5_old = -k_s_old * K_zz(ST(3)) * dx(1) / 2 * ((h_old(south) - h_old(point))/dz - 1);
gamma_6_old = -k_s_old * K_zz(ST(4)) * dx(2) / 2 * ((h_old(south) - h_old(point))/dz - 1);

% evaluate f function
f = phi(point) - phi_old(point) + dt/cell_volume * (theta * ...
            (gamma_1 + gamma_2 + gamma_3 + gamma_4 + gamma_5 + gamma_6) ...
            + (1 - theta) * ( ...
            gamma_1_old + gamma_2_old + gamma_3_old + gamma_4_old + gamma_5_old + gamma_6_old));
        
end
