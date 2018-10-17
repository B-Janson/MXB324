function [DIM]=VERIF_GRID_COORD()
% GRIDCOORD returns the dimension of the grid

% Width and height of aquifer
WIDTH = 500;
HEIGHT = 80;

% number of horizontal node points
n = 3;
% number of vertical node points
m = 11;

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

%% Create the distance matrix
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

%% Soil parameters (from PDF)
% Set all cells to be Alluvium
DIM.SP=[2.6,  0.91,   0.01,   0.33,   0.143, 1.51; ...   % Alluvium
        0.08, 0.0159, 0.106,  0.4686, 1.04,  1.3954; ... % 
        3.9,  1.17,   0.0286, 0.3658, 2.8,   2.239];     % 

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

DIM.r=r;
DIM.b=b;
% DIM.r = 1:n*m;

% DIM.NT = NT;

DIM.XZ = DIM.XZ(r,:);
DIM.NT = NT(r);
DIM.DELTA = DIM.DELTA(r, :);
DIM.VOL = DIM.VOL(r, :);
end