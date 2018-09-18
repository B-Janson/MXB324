function [PHI]=WCONT(DIM,h)
% n=DIM.n;
% m=DIM.m;
% XY=DIM.XY;

%Sres
%Ssat
%Cres
%Csat
Ares=0.01;
Asat=0.33;
PHI=zeros(1,length(h));
for i=1:length(h)
if h(i) >= 0
    PHI(i)=Asat;
else   
    PHI(i)=Ares+(Asat-Ares)*(1+(-2.8*h(i))^1.51)^(-(1-1/1.51));    
end
end

end