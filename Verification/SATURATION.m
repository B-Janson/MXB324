function [SAT] = SATURATION(DIM, h, i)

SAT = 0;
n = DIM.n_const(i, :);
m = 1 - 1./n;
alpha_const = DIM.alpha(i, :);

cell_volume = DIM.VOL(i, :);

if h(i) < 0
    for cell = 1:4
        SAT = SAT + cell_volume(cell) * (1 + (-alpha_const(cell) * h(i))^n(cell))^-m(cell);
    end
    SAT = SAT / cell_volume(5);
else
    SAT = 1;
end

end