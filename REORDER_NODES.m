function [DIM, h, S, phi, k, F] = REORDER_NODES(DIM, jacobian, h, S, phi, k, F)

r = symrcm(jacobian);
% r = 1:DIM.n*DIM.m;
b = bandwidth(jacobian(r, r));
bandwidth_lost = bandwidth(jacobian) - bandwidth(jacobian(r, r))

DIM.r = r;
DIM.b = b;

DIM.XZ = DIM.XZ(r, :);
DIM.NT = DIM.NT(r);
DIM.DELTA = DIM.DELTA(r, :);
DIM.VOL = DIM.VOL(r, :);
DIM.K_xx = DIM.K_xx(r, :);
DIM.K_zz = DIM.K_zz(r, :);
DIM.phi_res = DIM.phi_res(r, :);
DIM.phi_sat = DIM.phi_sat(r, :);
DIM.alpha = DIM.alpha(r, :);
DIM.n_const = DIM.n_const(r, :);

h = h(r);
S = S(r);
phi = phi(r);
k = k(r);
F = F(r);

end