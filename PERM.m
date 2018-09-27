function [k] = PERM(DIM, h, S, i)

k = 0;
n = DIM.n_const(i, :);
m = 1 - 1./n;
cell_volume = DIM.cell_volume(i, :);

if h(i) < 0
    for cell = 1:4
        k = k + cell_volume(cell) * (sqrt(S(i)) * (1 - (1 - S(i)^(1/m(cell)))^m(cell))^2);
    end
    k = k / sum(cell_volume);
else
    k = 1;
end

end