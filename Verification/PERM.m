function k = PERM(DIM, h, S, i)
% PERM  returns the permeability at the given node

% initialise k to 0
k = 0;
% get the soil types around this node
ST = DIM.ST(i, :);
% get each cell volume at this node
cell_volume = DIM.VOL(i, :);

if h(i) < 0
    % go through each cell and add up each permeability based on the soil
    % properties
    for cell = 1:4
        type = ST(cell);
        n = DIM.n_const(type);
        m = 1 - 1 / n;
        k = k + cell_volume(cell) * (sqrt(S(i)) * (1 - (1 - S(i)^(1/m))^m)^2);
    end
    % average k
    k = k / cell_volume(5);
else
    k = 1;
end

end