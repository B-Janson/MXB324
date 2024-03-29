function []=SOL_VIS(DIM, fig, colour_type, fig_title, vec)
% SOL_VIS  displays the solution of the given vector at whichever time step
%   draws to fig object, uses colour_type to determine colour scheme,
%   fig_title is the title the draw on the figure, and vec is the vector of
%   data points to visualise

% get grid information
n=DIM.n;
m=DIM.m;
Z=DIM.z;
X=DIM.x;
% revert grid to 'standard' view
vec(DIM.r) = vec;

VEC=zeros(m,n);
C=n*m+1;

for j = m:-1:1
    for k = n:-1:1
        C=C-1;
        VEC(j,k)=vec(C);
    end
end

% draw the surface plot and contour plot
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
