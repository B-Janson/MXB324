function [SAT] = SATURATION(DIM, h, i)
% SATURATION  returns the saturation at the given node

% initialise saturation to 0
SAT = 0;
% get the n value for this node
n = DIM.n_const(i, :);
% calculate m based on n
m = 1 - 1./n;
% get the alpha value for this node
alpha_const = DIM.alpha(i, :);
% get each cell volume at this node
cell_volume = DIM.VOL(i, :);

if h(i) < 0
    % go through each cell and add up each saturation based on the soil
    % properties
    for cell = 1:4
        SAT = SAT + cell_volume(cell) * (1 + (-alpha_const(cell) * h(i))^n(cell))^-m(cell);
    end
    % average saturation
    SAT = SAT / cell_volume(5);
else
    SAT = 1;
end

end