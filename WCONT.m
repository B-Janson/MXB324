function [PHI]=WCONT(DIM, h, S, i)
% Water content function

PHI = 0;
phi_res = DIM.phi_res(i, :);
phi_sat = DIM.phi_sat(i, :);
cell_volume = DIM.cell_volume(i, :);

if h(i) < 0
    for cell = 1:4
        PHI = PHI + cell_volume(cell) * (phi_res(cell) + S(i) * (phi_sat(cell) - phi_res(cell)));
    end
    PHI = PHI / sum(cell_volume);
else
    for cell = 1:4
        PHI = PHI + cell_volume(cell) * phi_sat(cell);
    end
    PHI = PHI / sum(cell_volume);
end

end