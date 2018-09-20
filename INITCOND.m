function [h, phi, k, S]=INITCOND(DIM)
n = DIM.n;
m = DIM.m;
XY = DIM.XY;
%%Create the initial vector

h   = zeros(n*m, 1);
phi = zeros(n*m, 1);
k   = zeros(n*m, 1);
S   = zeros(n*m, 1);

hbot = -5;
htop = -10;

for i = 1:n*m
    h(i) = hbot + (htop - hbot) * XY(i, 2)/80;
    S(i) = SATURATION(DIM, h, i);
    phi(i) = WCONT(DIM, h, S, i);
    k(i) = PERM(DIM, h, S, i);
end

h(DIM.r) = h;
S(DIM.r) = S;
phi(DIM.r) = phi;
k(DIM.r) = k;

end