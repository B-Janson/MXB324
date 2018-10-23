function [PUMPERS]=PUMPEVAL(DIM,RF,dt)
%Evaluate the pump rate at each node point
PUMPERS=DIM.S_P.*RF.*dt;

end