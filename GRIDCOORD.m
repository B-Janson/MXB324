function [DIM]=GRIDCOORD()
% GRIDCOORD returns the dimension of the grid

% Width and height of aquifer
WIDTH = 500;
HEIGHT = 80;

% number of horizontal node points
n = 251;
% number of vertical node points
m = 41;

num_nodes = n * m;

% Keep it uniform for now

% The discretisation we need in the x coordinate
DIM.x = linspace(0, WIDTH, n); 
DIM.n = n;

DIM.dx = zeros(1, n-1);
for i = 1:n-1
    DIM.dx(i) = DIM.x(i + 1) - DIM.x(i);
end

DELTAX = zeros(num_nodes, 2);
for j = 1:m
    DELTAX((j-1) * n + 1, :) = DIM.dx(1) / 2;
    DELTAX(j * n, :) = DIM.dx(end) / 2;

    for i = 2:n-1
        DELTAX((j-1) * n + i, 1) = DIM.dx(i - 1) / 2;
        DELTAX((j-1) * n + i, 2) = DIM.dx(i) / 2;
    end
end

%The Discretisation we need in the y
DIM.z = linspace(0, HEIGHT, m);
DIM.m = m;

DIM.dz = zeros(1, m-1);
for i = 1:m-1
    DIM.dz(i) = DIM.z(i + 1) - DIM.z(i);
end

DELTAZ = zeros(num_nodes, 2);
for j = 1:n
    DELTAZ((j-1) * m + 1, :) = DIM.dz(1) / 2;
    DELTAZ(j * m, :) = DIM.dz(end) / 2;

    for i = 2:m-1
        DELTAZ((j-1) * m + i, 1) = DIM.dz(i - 1) / 2;
        DELTAZ((j-1) * m + i, 2) = DIM.dz(i) / 2;
    end
end

DIM.cell_volume = zeros(n*m, 4);
for i = 1:m
    for j = 1:n
        index = (i-1) * n + j;
        DIM.cell_volume(index, 1) = DELTAX(index, 2) * DELTAZ(index, 2);
        DIM.cell_volume(index, 2) = DELTAX(index, 1) * DELTAZ(index, 2);
        DIM.cell_volume(index, 3) = DELTAX(index, 1) * DELTAZ(index, 1);
        DIM.cell_volume(index, 4) = DELTAX(index, 2) * DELTAZ(index, 1);
    end
end

%Create coordinate vector
[X, Z] = meshgrid(DIM.x, DIM.z);
X = X';
Z = Z';
XY = [X(:), Z(:)];
DIM.XY = XY;

% Set node point constants
DIM.K_xx = zeros(num_nodes, 4);
DIM.K_zz = zeros(num_nodes, 4);
DIM.phi_res = zeros(num_nodes, 4);
DIM.phi_sat = zeros(num_nodes, 4);
DIM.alpha = zeros(num_nodes, 4);
DIM.n_const = zeros(num_nodes, 4);

% Set all cells to be Alluvium
for i = 1:num_nodes
    DIM.K_xx(i, :)     = 2.6;
    DIM.K_zz(i, :)     = 0.91;
    DIM.phi_res(i, :)  = 0.01;
    DIM.phi_sat(i, :)  = 0.33;
    DIM.alpha(i, :)    = 1.43;
    DIM.n_const(i, :)  = 1.51;
end

% Assign a node type to each vertex,
% Just the baby problem for now
% Add in the river nodes later
NT = zeros(num_nodes, 1);

%Bottom L corner
NT(1) = 1;

%Bottom Row
NT(2:DIM.n-1) = 2;

%Bottom R corner
NT(n)=3;

for i = n+1:num_nodes-1
    if ((XY(i,1) == 0) && (XY(i,2) > 0) && (XY(i,2) < 30))
        %Left Sandstone Boundary
        NT(i)=4;
    elseif ( (XY(i,1) > 0) && (XY(i,1) < 500) && (XY(i,2) > 0) && (XY(i,2) < 30) ) || ...
            ( (XY(i,1) > 350) && (XY(i,1) < 500) && (XY(i,2) > 0) && (XY(i,2) < 80) )
        % Sandstone Interior
        % Pump at (450, 10)
        if ((XY(i,1) == 450) && (XY(i,2) == 10))
            NT(i) = 5;
        % Evapotranspiration zone (z >= 76)
        elseif XY(i,2) >= 76
            NT(i)=6;
        % Normal sandstone interior point
        else
            NT(i)=7;
        end
    elseif ((XY(i,1) == 500) && (XY(i,2) > 0) && (XY(i,2) < 80))
        % Right sandstone boundary
        if XY(i,2) < 76
            NT(i)=8;
        % Evapotranspiration zone
        else
            NT(i)=9;
        end
    elseif (XY(i,2) == 30)
        % Sandstone/Alluvium & Boundary
        if (XY(i,1) == 0)
            NT(i)=10;
        elseif ((XY(i,1) > 0) && (XY(i,1) < 50))
            % Sandstone/Alluvium Interface
            NT(i)=11;
        elseif (XY(i,1) == 50)
            % Sandstone/Alluvium/Confining interface
            NT(i)=12;
        elseif ((XY(i,1) > 50) && (XY(i,1) < 350))
            % Sandstone/Confining horizontal interface
            NT(i)=13;
        elseif (XY(i,1) == 350)
            % Sandstone/Confining corner
            NT(i)=14;
        end
    elseif ((XY(i,1) == 0) && (XY(i,2) > 30) && (XY(i,2) < 80))
        % Alluvium left boundary
        if XY(i,2) < 78
            NT(i)=15;
        else
            % Alluvium boundary evapotranspiration
            NT(i)=16;
        end
    elseif ( (XY(i,1) > 0) && (XY(i,1) < 50) && (XY(i,2) > 30) && (XY(i,2) < 80) ) || ...
            ( (XY(i,1) > 0) && (XY(i,1) < 350) && (XY(i,2) > 40) && (XY(i,2) < 80) )
        % Alluvium Interior
        if ((XY(i,1) == 100) && (XY(i,2) == 50))
            % Pump at (100, 50)
            NT(i)=17;
        elseif XY(i,2) >= 78
            % Evapotranspiration
            NT(i)=18;
        else
            NT(i)=19;
        end
    elseif ((XY(i,1) == 50) && (XY(i,2) > 30) && (XY(i,2) < 40))
        % Alluvium/Confining interface
        NT(i)=20;
    elseif ((XY(i,1) > 50) && (XY(i,1) < 350) && (XY(i,2) > 30) && (XY(i,2) < 40))
        % Confining interior
        NT(i)=21;
    elseif ((XY(i,1) == 350) && (XY(i,2) > 30) && (XY(i,2) < 40))
        % Confining/Sandstone vertical interface
        NT(i)=22;
    elseif ((XY(i,1) == 50) && (XY(i,2) == 40))
        % Confining/Alluvium corner
        NT(i)=23;
    elseif ((XY(i,1) > 50) && (XY(i,1) < 350) && (XY(i,2) == 40))
        % Confining/Alluvium horizontal interface
        NT(i)=24;
    elseif ((XY(i,1) == 350) && (XY(i,2) == 40))
        % Confining/Alluvium/Sandstone upper right corner
        NT(i)=25;
    elseif ((XY(i,1) == 350) && (XY(i,2) > 40) && (XY(i,2) < 80))
        % Alluvium/Sandstone interface
        if XY(i,2) < 76
            NT(i)=26;
        elseif XY(i,2) < 78
            % First half of evapo
            NT(i)=27;
        else
            % Second half of evapo
            NT(i)=28;
        end
    elseif 1 == 0 %((XY(i,1) == 0) && ((XY(i,2) == 60))
        %Bed/River
        NT(i)=28;
    elseif 1 == 0 %((XY(i,1) == 0) && ((XY(i,2) > 60))) && ((XY(i,2) <80))
        %River Boundary
        
        
    elseif ((XY(i,1) == 0) && (XY(i,2) == 80))
        % Top left node point River/Rain, L bedrock rain for baby problem
        NT(i)=29;
    elseif ((XY(i,1) > 0) && (XY(i,1) < 350) && (XY(i,2) == 80))
        % Alluvium/Rain boundary
        NT(i)=30;
    elseif ((XY(i,1) == 350) && (XY(i,2) == 80))
        % Alluvium/Sandstone/Rain Boundary
        NT(i)=31;
    elseif ((XY(i,1) > 350) && (XY(i,1) < 500) && (XY(i,2) == 80))
        % Sandstone/Rain Boundary
        NT(i)=32;
    end
end

NT(num_nodes)=33;

B=gallery('tridiag',num_nodes,1,1,1);
L=n;
U=num_nodes;
for i=1:n*(m-1)
    B(L,i)=1;
    B(U-n,U)=1;
    L=L+1;
    U=U-1;
end

r=symrcm(B);
b=bandwidth(B(r,r));

DIM.r=r;
DIM.b=b;

DIM.XY=DIM.XY(r,:);
DIM.NT=NT(r);
end