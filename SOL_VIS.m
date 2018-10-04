function []=SOL_VIS(DIM, fig, colour_type, fig_title, vec)

n=DIM.n;
m=DIM.m;
Z=DIM.z;
X=DIM.x;
vec(DIM.r) = vec;

VEC=zeros(m,n);
C=n*m+1;

for j=m:-1:1
    for k=1:n
        C=C-1;
        VEC(j,k)=vec(C);
    end
end

figure(fig)
hold on
surf(X,Z,VEC,'EdgeColor','none')
contour(X,Z,VEC,[0 0],'r','linewidth',0.5);
colormap(colour_type);
colorbar
drawnow
hold off

xlabel('X (m)')
ylabel('Z (m)')
title(fig_title)

end
