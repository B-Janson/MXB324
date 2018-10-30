function [phi_avg] = PHI_AVG(DIM, phi)
%PHI_AVG Calculate the numerical average of the water content

% Initialise to 0
phi_avg = 0;

% Loop through every node point
for i = 1:DIM.n*DIM.m
    % At each node point add the volume by the value
    phi_avg = phi_avg + DIM.VOL(i, 5) * phi(i);
end

% Take the weighted average over the whole domain
phi_avg = phi_avg / (DIM.WIDTH * DIM.HEIGHT);

end
