function [PHIBAR,Err]=MODVAL1(DIM,h,h0,t,r,rain_type)
x=DIM.x;
Lx=max(x);
y=DIM.y;
Ly=max(y);
PHI0=WCONT(DIM,h0);
PHI0=sum(PHI0);

if rain_type == 0
%For constant rainfall baby problem

PHIBAR=-(r/Lx)*t+PHI0;

elseif rain_type == 1
%For the sinusoidal function

PHIBAR= -(r/Lx)*t+ (365/(2*pi))*sin((2*pi)/365 * t)/(Lx*Ly) +PHI0;

elseif rain_type == 2
%For polynomial best fit of DALBY rainfall data

PHIBAR= 0; %get polynomial best fit to rainfall data over 12months    
    
end

Err=abs(PHIBAR-WCONT(DIM,h));


end