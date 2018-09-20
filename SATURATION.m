function [SAT] = SATURATION(DIM, h, i)

n = DIM.n_const(i);
m = 1 - 1/n;
alpha = DIM.alpha(i);

if h(i) < 0
    SAT = (1 + (-alpha * h(i))^n)^-m;
else
    SAT = 1;
end

end