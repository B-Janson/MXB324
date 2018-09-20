function [PHI]=WCONT(DIM, h, S, i)
% n=DIM.n;
% m=DIM.m;
% XY=DIM.XY;

%Sres
%Ssat
%Cres
%Csat
Ares=0.01;
Asat=0.33;
% PHI=zeros(1,length(h));

if h(i) < 0
    PHI = Ares + (Asat - Ares) * S(i);
else
    PHI = Asat;
end

end