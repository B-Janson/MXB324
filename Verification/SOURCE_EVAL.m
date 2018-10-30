function [PUMPS, EVAPS] = SOURCE_EVAL(DIM, PARAMS, PSI, RF)
% Evaluate the pump rate at each node point

% Pump vector is simple matrix multiplication
% DIM.S_P is a vector of the rate at each point
PUMPS = -PARAMS.PUMPING * DIM.S_P / DIM.WIDTH;

% Evapo vector is matrix multiplication of whether it should be active and
% the rate as given in DIM.S_E
ON = PSI ./ DIM.PSISAT_AVG > 0.5;
EVAPS = -ON .* DIM.S_E * RF / DIM.WIDTH;

end
