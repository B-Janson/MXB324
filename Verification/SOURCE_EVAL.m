function [PUMPERS,EVAPERS]=SOURCE_EVAL(DIM,RF,dt)
%Evaluate the pump rate at each node point
PUMPERS=DIM.S_P.*RF.*dt;
EVAPERS=DIM.S_E.*RF.*dt;
end