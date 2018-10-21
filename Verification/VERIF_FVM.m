function [F, S, phi, k] = VERIF_FVM(DIM, h, h_old, S_old, phi_old, k_old, t, PARAMS)

n = DIM.n;
m = DIM.m;
F = zeros(n*m, 1);
S = S_old;
phi = phi_old;
k = k_old;

for i = 1:n*m
    S(i) = SATURATION(DIM, h, i);
    phi(i) = WATER_CONTENT(DIM, h, S, i);
    k(i) = PERM(DIM, h, S, i);
end

DIM.r
for i = DIM.r
    DIM.NT(i)
    switch DIM.NT(i)
        case 1
            F(i) = V1(DIM, h, h_old, phi, phi_old, k, k_old, PARAMS);
        case 2
            F(i) = V2(DIM, i, h, h_old, phi, phi_old, k, k_old, PARAMS);
        case 3
            F(i) = V3(DIM, h, h_old, phi, phi_old, k, k_old, PARAMS);
        case 4
            F(i) = V4(DIM, i, h, h_old, phi, phi_old, k, k_old, PARAMS);
        case 5
            F(i) = V5(DIM, i, h, h_old, phi, phi_old, k, k_old, PARAMS);
        case 6
            F(i) = V6(DIM, i, h, h_old, phi, phi_old, k, k_old, PARAMS);
        case 7
            F(i) = V7(DIM, i, h, h_old, phi, phi_old, k, k_old, t, PARAMS);
        case 8
            F(i) = V8(DIM, i, h, h_old, phi, phi_old, k, k_old, t, PARAMS);
        case 9
            F(i) = V9(DIM, i, h, h_old, phi, phi_old, k, k_old, t, PARAMS);
        otherwise
            error('error')
            F(i) = V7(DIM, i, h, h_old, phi, phi_old, k, k_old, PARAMS);
    end
end


end