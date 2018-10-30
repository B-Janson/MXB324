function [F, S, phi, k, PUMPERS, EVAPERS] = VERIF_FVM(DIM, h, h_old, S_old, phi_old, k_old, PUMPERS_old, EVAPERS_old, RT, PARAMS)

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
        case 5
            F(i) = V5(DIM, DIM.r(i), h, h_old, phi, phi_old, k, k_old, PARAMS);
        case 6
            F(i) = V6(DIM, DIM.r(i), h, h_old, phi, phi_old, k, k_old, PARAMS);
        case 7
            F(i) = V7(DIM, h, h_old, phi, phi_old, k, k_old, RT, PARAMS);
        case 8
            F(i) = V8(DIM, DIM.r(i), h, h_old, phi, phi_old, k, k_old, RT, PARAMS);
        case 9
            F(i) = V9(DIM, h, h_old, phi, phi_old, k, k_old, RT, PARAMS);
        otherwise
            disp('Error! Wrong node type, assuming interior')
            F(i) = V7(DIM, i, h, h_old, phi, phi_old, k, k_old, PARAMS);
    end

end



%Deal with source terms

[PUMPERS,EVAPERS]=SOURCE_EVAL(DIM,RT,PARAMS.dt); 
sauce=PUMPERS+EVAPERS;
sauce_old=PUMPERS_old+EVAPERS_old;

F=F-(PARAMS.theta*sauce+(1-PARAMS.theta)*sauce_old);







end