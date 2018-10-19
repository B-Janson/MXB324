function [phi_avg] = PHI_AVG(DIM, phi)
%PHI_AVG Calculate the numerical average of the water content

phi_avg = 0;
for i = 1:DIM.n*DIM.m
    phi_avg = phi_avg + DIM.VOL(i, 5) * phi(i);
end

phi_avg = phi_avg / (DIM.WIDTH * DIM.HEIGHT);

end
