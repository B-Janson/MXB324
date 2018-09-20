function []=HEADVIS(DIM,PHI)

n=DIM.n;
m=DIM.m;
Z=DIM.z;
X=DIM.x;

HEAD=zeros(m,n);
C=n*m+1;

for j=m:-1:1
    for k=1:n
        C=C-1;
        HEAD(j,k)=PHI(C);
    end
end

figure()
hold on
surf(X,Z,HEAD,'EdgeColor','none')
contour(X,Z,HEAD,[0 0],'r','linewidth',0.5);
colorbar

xlabel('X (m)')
ylabel('Z (m)')
title('Water Content')