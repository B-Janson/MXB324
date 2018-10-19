function [SAT] = SATURATION(DIM, h, i)
% SATURATION  returns the saturation at the given node

% initialise saturation to 0
SAT = 0;
ST = DIM.ST(i, :);
% get each cell volume at this node
cell_volume = DIM.VOL(i, :);

if h(i) < 0
    % go through each cell and add up each saturation based on the soil
    % properties
    for cell = 1:4
        type = ST(cell);
        n = DIM.n_const(type);
        m = 1 - 1 / n;
        alpha_const = DIM.alpha(type);
        SAT = SAT + cell_volume(cell) * (1 + (-alpha_const * h(i))^n)^-m;
    end
    % average saturation
    SAT = SAT / cell_volume(5);
else
    SAT = 1;
end

end