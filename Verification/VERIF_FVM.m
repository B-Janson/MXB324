function [F, S, phi, k, PUMPS, EVAPOTS, river, total_rain] = VERIF_FVM(DIM, h, h_old, S_old, phi_old, k_old, t, PARAMS)
%VERIF_FVM function to calculate the finite volume method for each cell
%over the domain
% INPUTS
% DIM - struct of the dimension of the domain
% h - vector of pressure heads
% h_old - vector of pressure heads at the previous time step
% S_old - vector of saturation values at the previous time step
% phi_old - vector of water content values at the previous time step
% (should be called psi)
% k_old - vector of permeability values at the previous time step
% t - current time (days)
% PARAMS - struct of all problem parameters
% 
% OUTPUTS
% F - vector of function evaluations for each node/cell
% S - vector of saturation values for the given head/time
% phi - vectir of water content values for the given head/time
% k - vector of permeability values for the given head/time
% PUMPS - vector of amount of water has been lost to the pumps
% EVAPOTS - vector of amount of water has been lost to evaporation

% Get the dimensions of the domain
n = DIM.n;
m = DIM.m;
% Initialise the F output vector
F = zeros(n*m, 1);
S = S_old;
phi = phi_old;
k = k_old;
% Get the value of q dot n at the surface boundary
r_f = RAINFALL(PARAMS, t);
total_rain = r_f * PARAMS.dt / DIM.HEIGHT;
% initialise amount lost to water
river = 0;

% Calculate the new estimates for the closure conditions using this head
for i = 1:n*m
    S(i) = SATURATION(DIM, h, i);
    phi(i) = WATER_CONTENT(DIM, h, S, i);
    k(i) = PERM(DIM, h, S, i);
end

% Loop through each cell/node and return the evaluated internal functions
for i = 1:n*m
    switch DIM.NT(i)
        case 1
            F(i) = V1(DIM, h, h_old, phi, phi_old, k, k_old, r_f, PARAMS);
        case 2
            F(i) = V2(DIM, DIM.r(i), h, h_old, phi, phi_old, k, k_old, r_f, PARAMS);
        case 3
            F(i) = V3(DIM, h, h_old, phi, phi_old, k, k_old, r_f, PARAMS);
        case 4
            [F(i), node_river] = V4(DIM, DIM.r(i), h, h_old, phi, phi_old, k, k_old, r_f, PARAMS);
            river = river + node_river;
        case 5
            F(i) = V5(DIM, DIM.r(i), h, h_old, phi, phi_old, k, k_old, r_f, PARAMS);
        case 6
            F(i) = V6(DIM, DIM.r(i), h, h_old, phi, phi_old, k, k_old, r_f, PARAMS);
        case 7
            F(i) = V7(DIM, DIM.r(i), h, h_old, phi, phi_old, k, k_old, r_f, PARAMS);
        case 8
            F(i) = V8(DIM, DIM.r(i), h, h_old, phi, phi_old, k, k_old, r_f, PARAMS);
        case 9
            F(i) = V9(DIM, DIM.r(i), h, h_old, phi, phi_old, k, k_old, r_f, PARAMS);
        otherwise
            error('Error. Node point not set correctly at i = %d\n', i);
    end
end

river = river * PARAMS.dt / (DIM.WIDTH * DIM.HEIGHT);

% Source terms
[PUMPS, EVAPOTS] = SOURCE_EVAL(DIM, PARAMS, phi, r_f);
total_source = (PUMPS + EVAPOTS);

F = F - (PARAMS.theta * total_source + (1 - PARAMS.theta) * total_source);

end
