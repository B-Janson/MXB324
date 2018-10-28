function []=VISUALISE(DIM, h, psi, T, psi_true, psi_avg)
% SOL_VIS  displays the solution of the given vector at whichever time step
%   draws to fig object, uses colour_type to determine colour scheme,
%   fig_title is the title the draw on the figure, and vec is the vector of
%   data points to visualise

% get grid information
n = DIM.n;
m = DIM.m;
Z = DIM.z;
X = DIM.x;

% revert grid to 'standard' view
h(DIM.r) = h;
psi(DIM.r) = psi;

H = zeros(m, n);
PHI = zeros(m, n);
C = n * m + 1;

for j = m:-1:1
    for k = n:-1:1
        C = C - 1;
        H(j, k) = h(C);
        PHI(j, k) = psi(C);
    end
end

% draw the surface plot and contour plot
subplot(1, 3, 1)
hold on
surf(X, Z, H, 'EdgeColor','none')
contour(X, Z, H, [0 0], 'r', 'linewidth', 0.5);
colormap(gca, 'gray');
colorbar
drawnow
hold off
xlabel('X (m)')
ylabel('Z (m)')
title('Pressure Head (m)')

subplot(1, 3, 2)
hold on
surf(X, Z, PHI, 'EdgeColor','none')
colormap(gca, SAT_COLOUR);
colorbar
drawnow
hold off
xlabel('X (m)')
ylabel('Z (m)')
title('Water content')

subplot(1, 3, 3)
hold on
plot(T/365, psi_avg, 'b')
plot(T/365, psi_true, 'r')
legend('Simulation', 'Analytic', 'Location', 'southeast')
title('Water Content (Simulation vs Analytic)')
xlabel('Time (Years)')
ylabel('Water Content')
hold off
drawnow

end
