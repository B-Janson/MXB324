function [S]=SGEN(DIM,T,R)
n=DIM.n;
m=DIM.m;
XY=DIM.XY;

%% Create Source Vector

S=zeros(n*m,1);

%Low Pump
S(XY(1,:) == 450,XY(2,:) == 10)=; %mm/day
%High Pump
S(XY(1,:) == 100,XY(2,:) == 50)=; %mm/day


%Loop to determine evapotranspiration
%Change Start point of loop to improve efficency

for i=1:m*n
    if (XY(i,1) <= 350) && (XY(i,2) >= 78)
        S(i)=;
    elseif (XY(i,1) >= 350) && (XY(i,2) >= 76)
        S(i)=;
    end
end
end