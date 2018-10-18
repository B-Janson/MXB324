function [DIM]=GRIDCOORD()
% GRIDCOORD returns the dimension of the grid

% Width and height of aquifer
WIDTH = 500;
HEIGHT = 80;

% number of horizontal node points
n = 11;
% number of vertical node points
m = 9;

num_nodes = n * m;

% Keep it uniform for now

% Discretisation in x
DIM.x=linspace(0, WIDTH, n); 
DIM.n=n;

% Discretisation in z
DIM.z=linspace(0, HEIGHT, m);
DIM.m=m;

%Create coordinate vector
[X,Z]=meshgrid(DIM.x,DIM.z);
X=X';
Z=Z';
XZ=[X(:),Z(:)];
DIM.XZ=XZ;

% Create the distance matrix

% Create distance vectors dx=[L, R]
dx=zeros(DIM.n,2);
dx(1,2)=DIM.x(2)-DIM.x(1);
for i=2:DIM.n-1
    dx(i,1)=DIM.x(i)-DIM.x(i-1);
    dx(i,2)=DIM.x(i+1)-DIM.x(i);
end
dx(DIM.n,1)=DIM.x(DIM.n)-DIM.x(DIM.n-1);

% Create distance vectors dz=[D, U]
dz=zeros(DIM.m,2);
dz(1,2)=DIM.z(2)-DIM.z(1);
for i=2:DIM.m-1
    dz(i,1)=DIM.z(i)-DIM.z(i-1);
    dz(i,2)=DIM.z(i+1)-DIM.z(i);
end
dz(DIM.m,1)=DIM.z(DIM.m)-DIM.z(DIM.m-1);

% Create the big distances matrix
DELTA=zeros(num_nodes, 4);

c=0;
for i=0:n:n*(m-1)
    DELTA(i+1:i+n,1:2)=dx;
    c=c+1;
    for j=1:n
        DELTA(i+j,3:4)=dz(c,:);
    end
end

% Each row [L,R,D,U]
DIM.DELTA=DELTA;

% Each row contains the volumes
% [UR,UL,DL,DR,TV]
DIM.VOL = zeros(n*m, 5);
for i = 1:m*n
    DIM.VOL(i, 1) = DELTA(i, 2) * DELTA(i, 4) / 4; % UR
    DIM.VOL(i, 2) = DELTA(i, 1) * DELTA(i, 4) / 4; % UL
    DIM.VOL(i, 3) = DELTA(i, 1) * DELTA(i, 3) / 4; % DL
    DIM.VOL(i, 4) = DELTA(i, 2) * DELTA(i, 3) / 4; % DR
    DIM.VOL(i,5)=sum(DIM.VOL(i,1:4));
end

% Set node point constants
DIM.K_xx = zeros(num_nodes, 4);
DIM.K_zz = zeros(num_nodes, 4);
DIM.phi_res = zeros(num_nodes, 4);
DIM.phi_sat = zeros(num_nodes, 4);
DIM.alpha = zeros(num_nodes, 4);
DIM.n_const = zeros(num_nodes, 4);

% Set all cells to be Sandstone
for i = 1:num_nodes
    DIM.K_xx(i, :)     = 3.9;
    DIM.K_zz(i, :)     = 1.17;
    DIM.phi_res(i, :)  = 0.0286;
    DIM.phi_sat(i, :)  = 0.3658;
    DIM.alpha(i, :)    = 2.8;
    DIM.n_const(i, :)  = 2.239;
end

for i = 1:num_nodes
    x = XZ(i, 1);
    z = XZ(i, 2);
    if 0 <= x && x < 50 && z == 30
        % Alluvium top
        DIM.K_xx(i, 1:2)     = 2.6;
        DIM.K_zz(i, 1:2)     = 0.91;
        DIM.phi_res(i, 1:2)  = 0.01;
        DIM.phi_sat(i, 1:2)  = 0.33;
        DIM.alpha(i, 1:2)    = 1.43;
        DIM.n_const(i, 1:2)  = 1.51;
        
        % Sandstone bottom
        DIM.K_xx(i, 3:4)     = 3.9;
        DIM.K_zz(i, 3:4)     = 1.17;
        DIM.phi_res(i, 3:4)  = 0.0286;
        DIM.phi_sat(i, 3:4)  = 0.3658;
        DIM.alpha(i, 3:4)    = 2.8;
        DIM.n_const(i, 3:4)  = 2.239;
    elseif x == 50 && z == 30
        % Confining top right
        DIM.K_xx(i, 1)     = 0.08;
        DIM.K_zz(i, 1)     = 0.0159;
        DIM.phi_res(i, 1)  = 0.106;
        DIM.phi_sat(i, 1)  = 0.4686;
        DIM.alpha(i, 1)    = 1.04;
        DIM.n_const(i, 1)  = 1.3954;
        
        % Alluvium top left
        DIM.K_xx(i, 2)     = 2.6;
        DIM.K_zz(i, 2)     = 0.91;
        DIM.phi_res(i, 2)  = 0.01;
        DIM.phi_sat(i, 2)  = 0.33;
        DIM.alpha(i, 2)    = 1.43;
        DIM.n_const(i, 2)  = 1.51;
        
        % Sandstone bottom
        DIM.K_xx(i, 3:4)     = 3.9;
        DIM.K_zz(i, 3:4)     = 1.17;
        DIM.phi_res(i, 3:4)  = 0.0286;
        DIM.phi_sat(i, 3:4)  = 0.3658;
        DIM.alpha(i, 3:4)    = 2.8;
        DIM.n_const(i, 3:4)  = 2.239;
    elseif 50 < x && x < 350 && z == 30
        % Confining top
        DIM.K_xx(i, 1:2)     = 0.08;
        DIM.K_zz(i, 1:2)     = 0.0159;
        DIM.phi_res(i, 1:2)  = 0.106;
        DIM.phi_sat(i, 1:2)  = 0.4686;
        DIM.alpha(i, 1:2)    = 1.04;
        DIM.n_const(i, 1:2)  = 1.3954;
        
        % Sandstone bottom
        DIM.K_xx(i, 3:4)     = 3.9;
        DIM.K_zz(i, 3:4)     = 1.17;
        DIM.phi_res(i, 3:4)  = 0.0286;
        DIM.phi_sat(i, 3:4)  = 0.3658;
        DIM.alpha(i, 3:4)    = 2.8;
        DIM.n_const(i, 3:4)  = 2.239;
    elseif x == 350 && z == 30
        % Sandstone top right
        DIM.K_xx(i, 1)     = 3.9;
        DIM.K_zz(i, 1)     = 1.17;
        DIM.phi_res(i, 1)  = 0.0286;
        DIM.phi_sat(i, 1)  = 0.3658;
        DIM.alpha(i, 1)    = 2.8;
        DIM.n_const(i, 1)  = 2.239;
        
        % Confining top left
        DIM.K_xx(i, 2)     = 0.08;
        DIM.K_zz(i, 2)     = 0.0159;
        DIM.phi_res(i, 2)  = 0.106;
        DIM.phi_sat(i, 2)  = 0.4686;
        DIM.alpha(i, 2)    = 1.04;
        DIM.n_const(i, 2)  = 1.3954;
        
        % Sandstone bottom
        DIM.K_xx(i, 3:4)     = 3.9;
        DIM.K_zz(i, 3:4)     = 1.17;
        DIM.phi_res(i, 3:4)  = 0.0286;
        DIM.phi_sat(i, 3:4)  = 0.3658;
        DIM.alpha(i, 3:4)    = 2.8;
        DIM.n_const(i, 3:4)  = 2.239;
    elseif x == 350 && 30 < z && z < 40
        % Sandstone top right
        DIM.K_xx(i, 1)     = 3.9;
        DIM.K_zz(i, 1)     = 1.17;
        DIM.phi_res(i, 1)  = 0.0286;
        DIM.phi_sat(i, 1)  = 0.3658;
        DIM.alpha(i, 1)    = 2.8;
        DIM.n_const(i, 1)  = 2.239;
        
        % Confining left
        DIM.K_xx(i, 2:3)     = 0.08;
        DIM.K_zz(i, 2:3)     = 0.0159;
        DIM.phi_res(i, 2:3)  = 0.106;
        DIM.phi_sat(i, 2:3)  = 0.4686;
        DIM.alpha(i, 2:3)    = 1.04;
        DIM.n_const(i, 2:3)  = 1.3954;
        
        % Sandstone bottom right
        DIM.K_xx(i, 4)     = 3.9;
        DIM.K_zz(i, 4)     = 1.17;
        DIM.phi_res(i, 4)  = 0.0286;
        DIM.phi_sat(i, 4)  = 0.3658;
        DIM.alpha(i, 4)    = 2.8;
        DIM.n_const(i, 4)  = 2.239;
    elseif x == 350 && z == 40
        % Sandstone top right
        DIM.K_xx(i, 1)     = 3.9;
        DIM.K_zz(i, 1)     = 1.17;
        DIM.phi_res(i, 1)  = 0.0286;
        DIM.phi_sat(i, 1)  = 0.3658;
        DIM.alpha(i, 1)    = 2.8;
        DIM.n_const(i, 1)  = 2.239;
        
        % Alluvium top left
        DIM.K_xx(i, 2)     = 2.6;
        DIM.K_zz(i, 2)     = 0.91;
        DIM.phi_res(i, 2)  = 0.01;
        DIM.phi_sat(i, 2)  = 0.33;
        DIM.alpha(i, 2)    = 1.43;
        DIM.n_const(i, 2)  = 1.51;
        
        % Confining bottom left
        DIM.K_xx(i, 3)     = 0.08;
        DIM.K_zz(i, 3)     = 0.0159;
        DIM.phi_res(i, 3)  = 0.106;
        DIM.phi_sat(i, 3)  = 0.4686;
        DIM.alpha(i, 3)    = 1.04;
        DIM.n_const(i, 3)  = 1.3954;
        
        % Sandstone bottom right
        DIM.K_xx(i, 4)     = 3.9;
        DIM.K_zz(i, 4)     = 1.17;
        DIM.phi_res(i, 4)  = 0.0286;
        DIM.phi_sat(i, 4)  = 0.3658;
        DIM.alpha(i, 4)    = 2.8;
        DIM.n_const(i, 4)  = 2.239;
    elseif x == 350 && 40 < z && z < 80
        % Sandstone top right
        DIM.K_xx(i, 1)     = 3.9;
        DIM.K_zz(i, 1)     = 1.17;
        DIM.phi_res(i, 1)  = 0.0286;
        DIM.phi_sat(i, 1)  = 0.3658;
        DIM.alpha(i, 1)    = 2.8;
        DIM.n_const(i, 1)  = 2.239;
        
        % Alluvium left
        DIM.K_xx(i, 2:3)     = 2.6;
        DIM.K_zz(i, 2:3)     = 0.91;
        DIM.phi_res(i, 2:3)  = 0.01;
        DIM.phi_sat(i, 2:3)  = 0.33;
        DIM.alpha(i, 2:3)    = 1.43;
        DIM.n_const(i, 2:3)  = 1.51;
        
        % Sandstone bottom right
        DIM.K_xx(i, 4)     = 3.9;
        DIM.K_zz(i, 4)     = 1.17;
        DIM.phi_res(i, 4)  = 0.0286;
        DIM.phi_sat(i, 4)  = 0.3658;
        DIM.alpha(i, 4)    = 2.8;
        DIM.n_const(i, 4)  = 2.239;
    elseif x == 350 && z == 80
        % Alluvium left
        DIM.K_xx(i, 3)     = 2.6;
        DIM.K_zz(i, 3)     = 0.91;
        DIM.phi_res(i, 3)  = 0.01;
        DIM.phi_sat(i, 3)  = 0.33;
        DIM.alpha(i, 3)    = 1.43;
        DIM.n_const(i, 3)  = 1.51;
        
        % Sandstone bottom right
        DIM.K_xx(i, 4)     = 3.9;
        DIM.K_zz(i, 4)     = 1.17;
        DIM.phi_res(i, 4)  = 0.0286;
        DIM.phi_sat(i, 4)  = 0.3658;
        DIM.alpha(i, 4)    = 2.8;
        DIM.n_const(i, 4)  = 2.239;
    elseif x == 50 && 30 < z && z < 40
        % Confining top right
        DIM.K_xx(i, 1)     = 0.08;
        DIM.K_zz(i, 1)     = 0.0159;
        DIM.phi_res(i, 1)  = 0.106;
        DIM.phi_sat(i, 1)  = 0.4686;
        DIM.alpha(i, 1)    = 1.04;
        DIM.n_const(i, 1)  = 1.3954;
        
        % Alluvium left
        DIM.K_xx(i, 2:3)     = 2.6;
        DIM.K_zz(i, 2:3)     = 0.91;
        DIM.phi_res(i, 2:3)  = 0.01;
        DIM.phi_sat(i, 2:3)  = 0.33;
        DIM.alpha(i, 2:3)    = 1.43;
        DIM.n_const(i, 2:3)  = 1.51;
        
        % Confining bottom right
        DIM.K_xx(i, 4)     = 0.08;
        DIM.K_zz(i, 4)     = 0.0159;
        DIM.phi_res(i, 4)  = 0.106;
        DIM.phi_sat(i, 4)  = 0.4686;
        DIM.alpha(i, 4)    = 1.04;
        DIM.n_const(i, 4)  = 1.3954;
    elseif x == 50 && z == 40
        % Alluvium top and bottom left
        DIM.K_xx(i, 1:3)     = 2.6;
        DIM.K_zz(i, 1:3)     = 0.91;
        DIM.phi_res(i, 1:3)  = 0.01;
        DIM.phi_sat(i, 1:3)  = 0.33;
        DIM.alpha(i, 1:3)    = 1.43;
        DIM.n_const(i, 1:3)  = 1.51;
        
        % Confining bottom right
        DIM.K_xx(i, 4)     = 0.08;
        DIM.K_zz(i, 4)     = 0.0159;
        DIM.phi_res(i, 4)  = 0.106;
        DIM.phi_sat(i, 4)  = 0.4686;
        DIM.alpha(i, 4)    = 1.04;
        DIM.n_const(i, 4)  = 1.3954;
    elseif 50 < x && x < 350 && z == 40
        % Alluvium top
        DIM.K_xx(i, 1:2)     = 2.6;
        DIM.K_zz(i, 1:2)     = 0.91;
        DIM.phi_res(i, 1:2)  = 0.01;
        DIM.phi_sat(i, 1:2)  = 0.33;
        DIM.alpha(i, 1:2)    = 1.43;
        DIM.n_const(i, 1:2)  = 1.51;
        
        % Confining bottom
        DIM.K_xx(i, 3:4)     = 0.08;
        DIM.K_zz(i, 3:4)     = 0.0159;
        DIM.phi_res(i, 3:4)  = 0.106;
        DIM.phi_sat(i, 3:4)  = 0.4686;
        DIM.alpha(i, 3:4)    = 1.04;
        DIM.n_const(i, 3:4)  = 1.3954;
    elseif 0 <= x && x < 50 && 30 < z && z <= 80
        % Alluvium everywhere
        DIM.K_xx(i, :)     = 2.6;
        DIM.K_zz(i, :)     = 0.91;
        DIM.phi_res(i, :)  = 0.01;
        DIM.phi_sat(i, :)  = 0.33;
        DIM.alpha(i, :)    = 1.43;
        DIM.n_const(i, :)  = 1.51;
    elseif 50 <= x && x < 350 && 40 < z && z <= 80
        % Alluvium everywhere
        DIM.K_xx(i, :)     = 2.6;
        DIM.K_zz(i, :)     = 0.91;
        DIM.phi_res(i, :)  = 0.01;
        DIM.phi_sat(i, :)  = 0.33;
        DIM.alpha(i, :)    = 1.43;
        DIM.n_const(i, :)  = 1.51;
    elseif 50 < x && x < 350 && 30 < z && z < 40
        % Confining everywhere
        DIM.K_xx(i, :)     = 0.08;
        DIM.K_zz(i, :)     = 0.0159;
        DIM.phi_res(i, :)  = 0.106;
        DIM.phi_sat(i, :)  = 0.4686;
        DIM.alpha(i, :)    = 1.04;
        DIM.n_const(i, :)  = 1.3954;
    end     
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
    if XZ(i,1) == 0
        if XZ(i,2) == HEIGHT
            NT(i) = 7;
        else
            NT(i) = 4;
        end
    elseif XZ(i,1) == WIDTH
        NT(i) = 6;
    elseif XZ(i,2) == HEIGHT
        NT(i) = 8;
    else
        NT(i) = 5;
    end
end

NT(num_nodes)=9;

B=gallery('tridiag',num_nodes,1,1,1);
L=n;
U=num_nodes;
for i = 1:n*(m-1)
    B(L,i)=1;
    B(U-n,U)=1;
    L=L+1;
    U=U-1;
end

r=symrcm(B);
b=bandwidth(B(r,r));
Weightloss=2*(bandwidth(B)-bandwidth(B(r,r)))

r = 1:n*m;

DIM.r=r;
DIM.b=b;

DIM.XZ = DIM.XZ(r,:);
DIM.NT = NT(r);
DIM.DELTA = DIM.DELTA(r, :);
DIM.VOL = DIM.VOL(r, :);
DIM.K_xx = DIM.K_xx(r, :);
DIM.K_zz = DIM.K_zz(r, :);
DIM.phi_res = DIM.phi_res(r, :);
DIM.phi_sat = DIM.phi_sat(r, :);
DIM.alpha = DIM.alpha(r, :);
DIM.n_const = DIM.n_const(r, :);
end