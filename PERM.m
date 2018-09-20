function [k] = PERM(DIM, h, S, i)

n = DIM.n_const(i);
m = 1 - 1/n;

if h(i) < 0
    k = sqrt(S(i)) * (1 - (1 - S(i)^(1/m))^m)^2;
else
    k = 1;
end

end