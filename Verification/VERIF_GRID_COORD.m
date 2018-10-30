<<<<<<< HEAD
function [DIM]=VERIF_GRID_COORD(PARAMS, BC)
=======
function [DIM]=VERIF_GRID_COORD(LIN,PUMPS,EVAPOT)
>>>>>>> master
% GRIDCOORD returns the dimension of the grid

% Width and height of aquifer
WIDTH = 500;
HEIGHT = 80;

DIM.WIDTH = WIDTH;
DIM.HEIGHT = HEIGHT;

<<<<<<< HEAD
% number of horizontal node points
n = PARAMS.n;
% number of vertical node points
m = PARAMS.m;
=======
>>>>>>> master

num_nodes = n * m;

% Keep it uniform for now

% Discretisation in x
<<<<<<< HEAD
DIM.x = linspace(0, WIDTH, n);
DIM.n = n;
=======
DIM.x=[0,5,10:10:40,60:10:340,360:10:490,495,WIDTH];
%Make sure the boundaries are included

DIM.x(end+1:end+14)=[40,45,48,50,52,55,60,340,345,348,350,352,355,360];

DIM.x(end+1)=WIDTH;

%Tidy up
x=unique(DIM.x);
DIM.x=x;
n=length(DIM.x);
if LIN(1) == 1
    DIM.x=0:10:WIDTH;
    n=length(DIM.x);%If you want the linear version
end
DIM.n=n;
>>>>>>> master


% Discretisation in z
<<<<<<< HEAD
DIM.z = linspace(0, HEIGHT, m);
DIM.m = m;

% Create coordinate vector
[X,Z] = meshgrid(DIM.x,DIM.z);
X = X';
Z = Z';
XZ = [X(:), Z(:)];
DIM.XZ = XZ;
=======
DIM.z=[5:10:35,55:10:75];


DIM.z(end+1)=HEIGHT;
DIM.z(end+1)=HEIGHT-1;
DIM.z(end+1)=HEIGHT-2.5;
DIM.z(end+1:end+7)=[37.5,40,42.5,45,47.5,50,52.5];
DIM.z(end+1)=2.5;
DIM.z(end+1)=1;
DIM.z(end+1)=0;

%Tidy up
z=unique(DIM.z);
DIM.z=z;
m=length(DIM.z);
if LIN(2) == 1
    DIM.z=0:2:80;
    m=length(DIM.z);
    %If you want the Linear version
end
DIM.m=m;

num_o_nodes = n * m

%Create coordinate vector
[X,Z]=meshgrid(DIM.x,DIM.z);
X=X';
Z=Z';
XZ=[X(:),Z(:)];
DIM.XZ=XZ;
>>>>>>> master

% Create the distance matrix

% Create distance vectors dx=[L, R]
dx = zeros(DIM.n, 2);
dx(1,2) = DIM.x(2) - DIM.x(1);
for i = 2:DIM.n-1
    dx(i,1) = DIM.x(i) - DIM.x(i-1);
    dx(i,2) = DIM.x(i+1) - DIM.x(i);
end
dx(DIM.n,1) = DIM.x(DIM.n) - DIM.x(DIM.n-1);

% Create distance vectors dz=[D, U]
dz = zeros(DIM.m,2);
dz(1,2) = DIM.z(2) - DIM.z(1);
for i=2:DIM.m-1
    dz(i,1) = DIM.z(i) - DIM.z(i-1);
    dz(i,2) = DIM.z(i+1) - DIM.z(i);
end
dz(DIM.m,1) = DIM.z(DIM.m) - DIM.z(DIM.m-1);

% Create the big distances matrix
DELTA=zeros(num_o_nodes, 4);

c=0;
for i=0:n:n*(m-1)
    DELTA(i+1:i+n,1:2)=dx;
    c=c+1;
    for j=1:n
        DELTA(i+j,3:4)=dz(c,:);
    end
end

% Each row [L,R,D,U]
DIM.DELTA = DELTA;

% Each row contains the volumes
% [UR,UL,DL,DR,TV]
DIM.VOL = zeros(n*m, 5);
for i = 1:m*n
    DIM.VOL(i, 1) = DELTA(i, 2) * DELTA(i, 4) / 4; % UR
    DIM.VOL(i, 2) = DELTA(i, 1) * DELTA(i, 4) / 4; % UL
    DIM.VOL(i, 3) = DELTA(i, 1) * DELTA(i, 3) / 4; % DL
    DIM.VOL(i, 4) = DELTA(i, 2) * DELTA(i, 3) / 4; % DR
    DIM.VOL(i,5) = sum(DIM.VOL(i,1:4));
end

DIM.K_xx = [2.6 0.08 3.9];
DIM.K_zz = [0.91 0.0159 1.17];
DIM.phi_res = [0.01 0.106 0.0286];
DIM.phi_sat = [0.33 0.4686 0.3658];
DIM.alpha = [1.43 1.04 2.8];
DIM.n_const = [1.51 1.3954 2.239];

alluvium = 1;
confining = 2;
sandstone = 3;

% Set node point constants
DIM.ST = zeros(num_o_nodes, 4);

% Set all cells to be Sandstone
DIM.ST(:, :) = sandstone;

for i = 1:num_o_nodes
    x = XZ(i, 1);
    z = XZ(i, 2);
    if 0 <= x && x < 50 && z == 30
        % Alluvium top
        DIM.ST(i, 1:2) = alluvium;
        
        % Sandstone bottom
        DIM.ST(i, 3:4) = sandstone;
    elseif x == 50 && z == 30
        % Confining top right
        DIM.ST(i, 1) = confining;
        
        % Alluvium top left
        DIM.ST(i, 2) = alluvium;
        
        % Sandstone bottom
        DIM.ST(i, 3:4) = sandstone;
    elseif 50 < x && x < 350 && z == 30
        % Confining top
        DIM.ST(i, 1:2) = confining;
        
        % Sandstone bottom
        DIM.ST(i, 3:4) = sandstone;
    elseif x == 350 && z == 30
        % Sandstone top right
        DIM.ST(i, 1) = sandstone;
        
        % Confining top left
        DIM.ST(i, 2) = confining;
        
        % Sandstone bottom
        DIM.ST(i, 3:4) = sandstone;
    elseif x == 350 && 30 < z && z < 40
        % Sandstone top right
        DIM.ST(i, 1) = sandstone;
        
        % Confining left
        DIM.ST(i, 2:3) = confining;
        
        % Sandstone bottom right
        DIM.ST(i, 4) = sandstone;
    elseif x == 350 && z == 40
        % Sandstone top right
        DIM.ST(i, 1) = sandstone;
        
        % Alluvium top left
        DIM.ST(i, 2) = alluvium;
        
        % Confining bottom left
        DIM.ST(i, 3) = confining;
        
        % Sandstone bottom right
        DIM.ST(i, 4) = sandstone;
    elseif x == 350 && 40 < z && z < 80
        % Sandstone top right
        DIM.ST(i, 1) = sandstone;
        
        % Alluvium left
        DIM.ST(i, 2:3) = alluvium;
        
        % Sandstone bottom right
        DIM.ST(i, 4) = sandstone;
    elseif x == 350 && z == 80
        % Alluvium bottom left
        DIM.ST(i, 3) = alluvium;
        
        % Sandstone bottom right
        DIM.ST(i, 4) = sandstone;
    elseif x == 50 && 30 < z && z < 40
        % Confining top right
        DIM.ST(i, 1) = confining;
        
        % Alluvium left
        DIM.ST(i, 2:3) = alluvium;
        
        % Confining bottom right
        DIM.ST(i, 4) = confining;
    elseif x == 50 && z == 40
        % Alluvium top and bottom left
        DIM.ST(i, 1:3) = alluvium;
        
        % Confining bottom right
        DIM.ST(i, 4) = confining;
    elseif 50 < x && x < 350 && z == 40
        % Alluvium top
        DIM.ST(i, 1:2) = alluvium;
        
        % Confining bottom
        DIM.ST(i, 3:4) = confining;
    elseif 0 <= x && x < 50 && 30 < z && z <= 80
        % Alluvium everywhere
        DIM.ST(i, :) = alluvium;
    elseif 50 <= x && x < 350 && 40 < z && z <= 80
        % Alluvium everywhere
        DIM.ST(i, :) = alluvium;
    elseif 50 < x && x < 350 && 30 < z && z < 40
        % Confining everywhere
        DIM.ST(i, :) = confining;
    end
end

% Assign a node type to each vertex,
% Just the baby problem for now
% Add in the river nodes later
<<<<<<< HEAD
NT = zeros(num_nodes, 1);

for i = 1:m*n
    if XZ(i, 2) == 0
        % Bottom edge
        NT(i) = 2;
    elseif XZ(i, 2) == HEIGHT
        % Top edge
        NT(i) = 8;
    elseif XZ(i, 1) == 0
        % Left edge
        NT(i) = 4;
    elseif XZ(i, 1) == WIDTH
        % Right edge
        NT(i) = 6;
=======
NT = zeros(num_o_nodes, 1);


for i=1:m*n
    if XZ(i,2) == 0
        %Bottom edge
        NT(i)=2;
    elseif XZ(i,2) == 80
        %Top edge
        NT(i)=8;
    elseif XZ(i,1) == 0
        %Left edge
        NT(i)=4;
    elseif XZ(i,1) == 500
        %Right edge
        NT(i)=6;
>>>>>>> master
    else
        % Interior
        NT(i) = 5;
    end
end
<<<<<<< HEAD

% Bottom Left
NT(1) = 1;
% Bottom Right
NT(n) = 3;
% Top Left
NT(n * (m - 1) + 1) = 7;
% Top Right
NT(n*m) = 9;

NODE_BC = zeros(num_nodes, 2);

NODE_BC(1, 1) = BC.left(1);
NODE_BC(1, 2) = BC.bottom(1);

NODE_BC(n, 1) = BC.right(1);
NODE_BC(n, 2) = BC.bottom(end);

NODE_BC(end - n + 1, 1) = 0; % rainfall
NODE_BC(end - n + 1, 2) = BC.left(end);

NODE_BC(end, 1) = BC.right(end);
NODE_BC(end, 2) = 0; % rainfall

for i = 2:n-1
    idx = 2 * (i - 1);
    NODE_BC(i, :) = BC.bottom(idx:idx+1);
end

for i = 2:m-1
    index = (i - 1) * n + 1;
    idx = 2 * (i - 1);
    NODE_BC(index, :) = BC.left(idx:idx+1);
    
    index = i * n;
    idx = 2 * (i - 1);
    NODE_BC(index, :) = BC.right(idx:idx+1);
end

% Discover layout of the jacobian
off_diag = ones(1, num_nodes - 1);
off_diag(n:n:end) = 0;
far_band = ones(1, num_nodes - n);

B = diag(ones(1, num_nodes)) + diag(off_diag, 1) + diag(off_diag, -1) ... 
    + diag(far_band, n) + diag(far_band, -n);


% RCM reorder the jacobian
r = symrcm(B);
Weightloss=2*(bandwidth(B)-bandwidth(B(r,r)));
b = 2 * bandwidth(B(r,r)) + 1;

DIM.r = r;
DIM.b = b;

% Reorder Everything
DIM.XZ = DIM.XZ(r, :);
=======
%Bottom Left
NT(1)=1;
%Bottom Right
NT(n)=3;
%Top Left
NT(n*(m-1)+1)=7;
%Top Right
NT(end)=9;

%Discover layout of the jacobian
B=gallery('tridiag',num_o_nodes,1,1,1);
L=n;
U=num_o_nodes;
for i = 1:n*(m-1)
    B(L,i)=1;
    B(U-n,U)=1;
    L=L+1;
    U=U-1;
end

%% Pump Terms
%Each pump takes up a row of PUMPS=[x,z,fraction of annual rainfall]


S_P=zeroes(num_o_nodes,1);
[LP,~]=size(PUMPS);
j=1;
while j <= LP
    for i=1:num_o_nodes        
        if (DIM.XZ(1,i) == PUMPS(j,1)) &&(DIM.XZ(1,i) == PUMPS(j,1))
            S_P(i)=PUMPS(3,j);
            j=j+1;
            break
        end
    end
end


%% Evapotranspiration Terms 
%Defines a region of evapotranspiration as such [L,R,Depth,fraction of annual rainfall] 

S_E=zeroes(num_o_nodes,1);
[LE,~]=size(EVAPOT);
j=1;

while j <= LE
    DEPTH=EVAPOT(j,3);
    for i=1:num_o_nodes        
        %Check if point is in jth evapotranspiration zone
        if ((DIM.XZ(1,i) <= EVAPOT(j,2)) &&(DIM.XZ(1,i) >= EVAPOT(j,1))) %Check X-direction
            if (DIM.XZ(2,i) >= HEIGHT-DEPTH) %Check Z-direction
            S_E(i)=-EVAPOT(j,4)*(DIM.z(i,2)-WIDTH+DEPTH)^2/(DEPTH)^2;%Q(z) in the evapotranspiration zone           
            end
        end
    end 
    j=j+1;
end



%Check how much water is beeing sucked out
if (sum(EVAPOT(:,4)+sum(PUMPS(:,3))) >= 1)
    disp('Caution! Too much SUCC')
end

%% Reordering

%RCM reorder the jacobian
r=symrcm(B);
b=bandwidth(B(r,r));
Weightloss=2*(bandwidth(B)-bandwidth(B(r,r)))
DIM.r=r;
DIM.b=b;

%Reorder Everything
DIM.S_P=S_P(r);
DIM.S_E=S_E(r);
DIM.PUMPS=PUMPS(r,:);
DIM.EVAPO=EVAPOT(r,:);
DIM.XZ = DIM.XZ(r,:);
>>>>>>> master
DIM.NT = NT(r);
DIM.ST = DIM.ST(r, :);
DIM.DELTA = DIM.DELTA(r, :);
DIM.VOL = DIM.VOL(r, :);
DIM.BC = NODE_BC(r, :);
end