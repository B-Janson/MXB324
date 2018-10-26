function [F, S, phi, k] = VERIF_FVM(DIM, h, h_old, S_old, phi_old, k_old, t, PARAMS)

n = DIM.n;
m = DIM.m;
F = zeros(n*m, 1);
S = S_old;
phi = phi_old;
k = k_old;
r_f = RAINFALL(PARAMS, t);

for i = 1:n*m
    S(i) = SATURATION(DIM, h, i);
    phi(i) = WATER_CONTENT(DIM, h, S, i);
    k(i) = PERM(DIM, h, S, i);
end

for i = 1:n*m
    switch DIM.NT(i)
        case 1
            F(i) = V1(DIM, h, h_old, phi, phi_old, k, k_old, r_f, PARAMS);
        case 2
            F(i) = V2(DIM, DIM.r(i), h, h_old, phi, phi_old, k, k_old, r_f, PARAMS);
        case 3
            F(i) = V3(DIM, h, h_old, phi, phi_old, k, k_old, r_f, PARAMS);
        case 4
            F(i) = V4(DIM, DIM.r(i), h, h_old, phi, phi_old, k, k_old, r_f, PARAMS);
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

% Source terms
[PUMPS, EVAPOTS] = SOURCE_EVAL(DIM, PARAMS, phi, r_f);
total_source = PUMPS + EVAPOTS;

F = F - PARAMS.dt * (PARAMS.theta * total_source + (1 - PARAMS.theta) * total_source);

end
