function f = V30(DIM, i, h, h_old, phi, phi_old, k, k_old, t, PARAMS)
% Sandstone/Rain

n = DIM.n;
D=DIM.DELTA(i,:);
K_xx = DIM.SP(DIM.ST(i,2),2);
K_zz = DIM.SP(DIM.ST(i,1),2);
VOL = DIM.VOL(i,5);

active_flow = 1;
active_flow_old = 1;
dt = PARAMS.dt;
theta = PARAMS.theta;
r_f = PARAMS.r_f;
point = (DIM.r == i);
east = (DIM.r == i+1);
west = (DIM.r == i-1);
south = (DIM.r == i-n);

% if h(point) >= 0
%    active_flow = 0;
% end

% if h_old(point) >= 0
%    active_flow_old = 0; 
% end



k_s = (k(point) + k(south)) / 2;
k_s_old = (k_old(point) + k_old(south)) / 2;

gamma_1 = -active_flow * r_f * D(1) / 2;
gamma_2 = k_s * K_zz * ((h(south) - h(point)) / D(3) + 1) * D(1) / 2;

gamma_1_old = -active_flow_old * r_f * D(2) / 2;
gamma_2_old = k_s_old * K_zz * ((h_old(south) - h_old(point)) / D(3) + 1) * D(1) / 2;

f = phi(point) - phi_old(point) + dt/VOL * (theta * ( ...
            gamma_1 + gamma_2) ...
            + (1 - theta) * ( ...
            gamma_1_old + gamma_2_old));
        
end
