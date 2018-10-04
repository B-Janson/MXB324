function F = FVM_TEST(DIM, h, h_old, S_old, phi_old, k_old, t, PARAMS)

n = DIM.n;
m = DIM.m;
F = zeros(n*m, 1);
% h_old(DIM.r) = h_old;
% h(DIM.r) = h;
% S_old(DIM.r) = S_old;
% phi_old(DIM.r) = phi_old;
% k_old(DIM.r) = k_old;
S = S_old;
phi = phi_old;
k = k_old;

r_f = 0.00171;

for i = 1:n*m
    S(i) = SATURATION(DIM, h, i);
    phi(i) = WATER_CONTENT(DIM, h, S, i);
    k(i) = PERM(DIM, h, S, i);
end

for i = 1:n*m
    switch DIM.NT(i)
        case 1
            F(i) = V1(DIM, h, h_old, phi, phi_old, k, k_old, PARAMS);
        case 2
            F(i) = V2(DIM, DIM.r(i), h, h_old, phi, phi_old, k, k_old, PARAMS);
        case 3
            F(i) = V3(DIM, h, h_old, phi, phi_old, k, k_old, PARAMS);
        case 4
            F(i) = V4(DIM, DIM.r(i), h, h_old, phi, phi_old, k, k_old, PARAMS);
        case 7
            F(i) = V7(DIM, DIM.r(i), h, h_old, phi, phi_old, k, k_old, PARAMS);
        case 8
            F(i) = V8(DIM, DIM.r(i), h, h_old, phi, phi_old, k, k_old, PARAMS);
        case 15
            F(i) = V15(DIM, DIM.r(i), h, h_old, phi, phi_old, k, k_old, PARAMS);
        case 19
            F(i) = V19(DIM, DIM.r(i), h, h_old, phi, phi_old, k, k_old, PARAMS);
        case 24
            F(i) = V24(DIM, DIM.r(i), h, h_old, phi, phi_old, k, k_old, PARAMS);
        case 29
            F(i) = V29(DIM, DIM.r(i), h, h_old, phi, phi_old, k, k_old, t, PARAMS);
        case 30
            F(i) = V30(DIM, DIM.r(i), h, h_old, phi, phi_old, k, k_old, t, PARAMS);
        case 31
            F(i) = V31(DIM, DIM.r(i), h, h_old, phi, phi_old, k, k_old, t, PARAMS);
        case 32
            F(i) = V32(DIM, DIM.r(i), h, h_old, phi, phi_old, k, k_old, t, PARAMS);
        case 33
            F(i) = V33(DIM, DIM.r(i), h, h_old, phi, phi_old, k, k_old, t, PARAMS);
        otherwise
            F(i) = V7(DIM, DIM.r(i), h, h_old, phi, phi_old, k, k_old, PARAMS);
    end
end


end