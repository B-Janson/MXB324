function [PUMPS, EVAPS] = SOURCE_EVAL(DIM, PARAMS, PSI, RF)
% Evaluate the pump rate at each node point

PUMPS = -PARAMS.PUMPING * DIM.S_P / DIM.WIDTH;

ON = PSI ./ DIM.PSISAT_AVG > 0.5;
EVAPS = -ON .* DIM.S_E * RF / DIM.WIDTH;

end
