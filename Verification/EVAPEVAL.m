function [EVAPOTRANSPIRATION]=EVAPEVAL(DIM,RF,dt)
%Evaluate the pump rate at each node point

for i=1:DIM.



EVAPOTRANSPIRATION=-DIM.S_E.*RF.*dt;

end