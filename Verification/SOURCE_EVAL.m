function [PUMPERS,EVAPERS]=SOURCE_EVAL(DIM,PARAMS,RF,dt)
%Evaluate the pump rate at each node point
PUMPERS=PARAMS.PUMPS.*DIM.S_P.*dt/DIM.WIDTH;

ON=(PHI > DIM.PHISAT); 
EVAPERS=ON.*DIM.S_E.*RF.*dt/DIM.WIDTH;

end