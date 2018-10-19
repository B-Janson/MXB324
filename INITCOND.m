function [h, S, phi, k]=INITCOND(DIM)
% Function to set up all variables according the initial state of the
% aquifer

n = DIM.n;
m = DIM.m;
XZ = DIM.XZ;

% Given initial constants
hbot = -5;
htop = -10;

% Create the initial vectors
h   = zeros(n*m, 1);
phi = zeros(n*m, 1);
k   = zeros(n*m, 1);
S   = zeros(n*m, 1);

for i = 1:n*m
    h(i) = hbot + (htop - hbot) * XZ(i, 2)/80;
    S(i) = SATURATION(DIM, h, i);
    phi(i) = WATER_CONTENT(DIM, h, S, i);
    k(i) = PERM(DIM, h, S, i);
end

% h(DIM.r) = h;
% S(DIM.r) = S;
% phi(DIM.r) = phi;
% k(DIM.r) = k;

end