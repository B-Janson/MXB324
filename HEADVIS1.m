function []=HEADVIS1(DIM,H)

n=DIM.n;
m=DIM.m;
Y=DIM.y;
X=DIM.x;

HEAD=zeros(m,n);
C=n*m+1;

for j=m:-1:1
    for k=1:n
        C=C-1;
        HEAD(j,k)=H(C);
    end
end


figure()
hold on
surf(X,Y,HEAD,'EdgeColor','none')
contour(X,Y,HEAD,[0 0],'r','linewidth',0.5);
colormap('gray')
colorbar

xlabel('X (m)')
ylabel('Z (m)')
title('Pressure Head (m)')



