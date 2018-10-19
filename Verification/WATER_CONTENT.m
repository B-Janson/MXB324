function [PHI]=WATER_CONTENT(DIM, h, S, i)
% Water content function

PHI = 0;
ST = DIM.ST(i, :);
cell_volume = DIM.VOL(i, :);

if h(i) < 0
    for cell = 1:4
        type = ST(cell);
        res = DIM.phi_res(type);
        sat = DIM.phi_sat(type);
        PHI = PHI + cell_volume(cell) * (res + S(i) * (sat - res));
    end
else
    for cell = 1:4
        type = ST(cell);
        sat = DIM.phi_sat(type);
        PHI = PHI + cell_volume(cell) * sat;
    end
end

PHI = PHI / cell_volume(5);

end