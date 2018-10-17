function [PHI]=WATER_CONTENT(DIM, h, S, i)
% Water content function

PHI = 0;
phi_res = DIM.phi_res(i, :);
phi_sat = DIM.phi_sat(i, :);
cell_volume = DIM.VOL(i, :);

if h(i) < 0
    for cell = 1:4
        PHI = PHI + cell_volume(cell) * (phi_res(cell) + S(i) * (phi_sat(cell) - phi_res(cell)));
    end
else
    for cell = 1:4
        PHI = PHI + cell_volume(cell) * phi_sat(cell);
    end
end

PHI = PHI / cell_volume(5);

end