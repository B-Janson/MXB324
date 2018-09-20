function [U0]=INITCOND(DIM)
n=DIM.n;
m=DIM.m;
XY=DIM.XY;
%%Create the initial vector

U0=zeros(n*m,1);
hbot=1;
htop=-10;



for i=1:n*m
    U0(i)=hbot+(htop-hbot)*XY(i,2)/80;
end
U0(DIM.r)=U0;
end