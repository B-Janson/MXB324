function [DIM] = VERIF_GRID_COORD(PARAMS)
% GRIDCOORD returns the dimension of the grid

% Width and height of aquifer
WIDTH = 500;
HEIGHT = 80;

DIM.WIDTH = WIDTH;
DIM.HEIGHT = HEIGHT;
PARAMS.uniform
if PARAMS.uniform == 1
    
    % Discretisation in x
    DIM.x = 0:25:WIDTH; %course grid
    
    % Discretisation in z
    DIM.z = 0:2:HEIGHT;
    
    n = length(DIM.x);
    m = length(DIM.z);
else
    % Discretisation in x
    DIM.x = [0,5,10:15:40,45,50,55,60:15:90,95,...
        100,105,linspace(110,340,9),345,350,355,...
        360:20:440,445,460:15:490,495,500];
    n = length(DIM.x)
    
    % Discretisation in z
    DIM.z = [0,1,2,4,8,12,16,20,24,28,29,30,31,32,34,36,...
        38,39,40,41,42,44,48,52,56,58,60,62,64,65,66,70,74,78,79,80];
    m = length(DIM.z)
end

% Ensure that all corner points are included
if sum(ismember(DIM.x, 0) | ismember(DIM.x, 500)) + sum(ismember(DIM.z, 0) | ismember(DIM.z, 80)) ~= 4
    error('Grid must contain (0, 0), (500, 0), (0, 80) & (500, 80)');
end

num_nodes = n * m;

DIM.n = n;
DIM.m = m;

% Create coordinate vector
[X,Z] = meshgrid(DIM.x, DIM.z);
X = X';
Z = Z';
XZ = [X(:), Z(:)];

% Create the distance matrix

% Create distance vectors dx = [L, R]
dx = zeros(n, 2);
dx(1,2) = DIM.x(2) - DIM.x(1);
for i = 2:n-1
    dx(i,1) = DIM.x(i) - DIM.x(i-1);
    dx(i,2) = DIM.x(i+1) - DIM.x(i);
end
dx(n, 1) = DIM.x(n) - DIM.x(n-1);

% Create distance vectors dz = [D, U]
dz = zeros(m, 2);
dz(1,2) = DIM.z(2) - DIM.z(1);
for i = 2:m-1
    dz(i,1) = DIM.z(i) - DIM.z(i-1);
    dz(i,2) = DIM.z(i+1) - DIM.z(i);
end
dz(m, 1) = DIM.z(m) - DIM.z(m-1);

% Create the big distances matrix
DELTA = zeros(num_nodes, 4);

c = 0;
for i = 0:n:n*(m-1)
    DELTA(i+1:i+n,1:2) = dx;
    c = c+1;
    for j = 1:n
        DELTA(i+j,3:4) = dz(c, :);
    end
end

% Each row contains the volumes
% [UR, UL, DL, DR, Total]
VOL = zeros(n*m, 5);
for i = 1:m*n
    VOL(i, 1) = DELTA(i, 2) * DELTA(i, 4) / 4; % UR
    VOL(i, 2) = DELTA(i, 1) * DELTA(i, 4) / 4; % UL
    VOL(i, 3) = DELTA(i, 1) * DELTA(i, 3) / 4; % DL
    VOL(i, 4) = DELTA(i, 2) * DELTA(i, 3) / 4; % DR
    VOL(i,5) = sum(VOL(i,1:4));
end

% All soil properties [Alluvium, Confining, Sandstone]
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
ST = zeros(num_nodes, 4);

% Set all cells to be Sandstone
ST(:, :) = sandstone;

for i = 1:num_nodes
    x = XZ(i, 1);
    z = XZ(i, 2);
    if 0 <= x && x < 50 && z == 30
        % Alluvium top
        ST(i, 1:2) = alluvium;
        
        % Sandstone bottom
        ST(i, 3:4) = sandstone;
    elseif x == 50 && z == 30
        % Confining top right
        ST(i, 1) = confining;
        
        % Alluvium top left
        ST(i, 2) = alluvium;
        
        % Sandstone bottom
        ST(i, 3:4) = sandstone;
    elseif 50 < x && x < 350 && z == 30
        % Confining top
        ST(i, 1:2) = confining;
        
        % Sandstone bottom
        ST(i, 3:4) = sandstone;
    elseif x == 350 && z == 30
        % Sandstone top right
        ST(i, 1) = sandstone;
        
        % Confining top left
        ST(i, 2) = confining;
        
        % Sandstone bottom
        ST(i, 3:4) = sandstone;
    elseif x == 350 && 30 < z && z < 40
        % Sandstone top right
        ST(i, 1) = sandstone;
        
        % Confining left
        ST(i, 2:3) = confining;
        
        % Sandstone bottom right
        ST(i, 4) = sandstone;
    elseif x == 350 && z == 40
        % Sandstone top right
        ST(i, 1) = sandstone;
        
        % Alluvium top left
        ST(i, 2) = alluvium;
        
        % Confining bottom left
        ST(i, 3) = confining;
        
        % Sandstone bottom right
        ST(i, 4) = sandstone;
    elseif x == 350 && 40 < z && z < 80
        % Sandstone top right
        ST(i, 1) = sandstone;
        
        % Alluvium left
        ST(i, 2:3) = alluvium;
        
        % Sandstone bottom right
        ST(i, 4) = sandstone;
    elseif x == 350 && z == 80
        % Alluvium bottom left
        ST(i, 3) = alluvium;
        
        % Sandstone bottom right
        ST(i, 4) = sandstone;
    elseif x == 50 && 30 < z && z < 40
        % Confining top right
        ST(i, 1) = confining;
        
        % Alluvium left
        ST(i, 2:3) = alluvium;
        
        % Confining bottom right
        ST(i, 4) = confining;
    elseif x == 50 && z == 40
        % Alluvium top and bottom left
        ST(i, 1:3) = alluvium;
        
        % Confining bottom right
        ST(i, 4) = confining;
    elseif 50 < x && x < 350 && z == 40
        % Alluvium top
        ST(i, 1:2) = alluvium;
        
        % Confining bottom
        ST(i, 3:4) = confining;
    elseif 0 <= x && x < 50 && 30 < z && z <= 80
        % Alluvium everywhere
        ST(i, :) = alluvium;
    elseif 50 <= x && x < 350 && 40 < z && z <= 80
        % Alluvium everywhere
        ST(i, :) = alluvium;
    elseif 50 < x && x < 350 && 30 < z && z < 40
        % Confining everywhere
        ST(i, :) = confining;
    end
end

% Assign a node type to each vertex,
NT = zeros(num_nodes, 1);

for i = 1:num_nodes
    x = XZ(i, 1);
    z = XZ(i, 2);
    if z == 0
        % Bottom edge
        NT(i) = 2;
    elseif z == HEIGHT
        % Top edge
        NT(i) = 8;
    elseif x == 0
        % Left edge
        NT(i) = 4;
    elseif x == WIDTH
        % Right edge
        NT(i) = 6;
    else
        % Interior
        NT(i) = 5;
    end
end

% Bottom Left
NT(1) = 1;
% Bottom Right
NT(n) = 3;
% Top Left
NT(n * (m - 1) + 1) = 7;
% Top Right
NT(n*m) = 9;

% Pump terms
% Each pump takes up a row of PUMPS=[x, z, rate]

S_P = zeros(num_nodes, 1);
[num_pumps, ~] = size(PARAMS.PUMPS);
for pump = 1:num_pumps
    for i = 1:num_nodes
        if XZ(i, 1) == PARAMS.PUMPS(pump, 1) && XZ(i, 2) == PARAMS.PUMPS(pump, 2)
            S_P(i) = PARAMS.PUMPS(pump, 3);
            break;
        end
    end
end

% Evapotranspiration
% Defines a region of evapotranspiration as such [L, R, Depth, fraction]

S_E = zeros(num_nodes, 1);
[num_zones, ~] = size(PARAMS.EVAPOT);

for zone = 1:num_zones
    current_zone = PARAMS.EVAPOT(zone, :);
    DEPTH = current_zone(3);
    for i = 1:num_nodes
        x = XZ(i, 1);
        z = XZ(i, 2);
        % Check if point is in this zone
        if current_zone(1) <= x && x <= current_zone(2)
            if z > HEIGHT - DEPTH
                S_E(1) = -current_zone(4) * (z - HEIGHT + DEPTH)^2 / (DEPTH)^2;
            end
        end
    end
end

PSISAT_AVG = zeros(num_nodes, 1);

for i = 1:num_nodes
    for j = 1:4
        PSISAT_AVG(i) = PSISAT_AVG(i) + VOL(i, j) * DIM.phi_sat(ST(i, j));
    end
end

PSISAT_AVG = PSISAT_AVG ./ VOL(:, 5);

% Discover layout of the jacobian
off_diag = ones(1, num_nodes - 1);
off_diag(n:n:end) = 0;
far_band = ones(1, num_nodes - n);

B = diag(ones(1, num_nodes)) + diag(off_diag, 1) + diag(off_diag, -1) ...
    + diag(far_band, n) + diag(far_band, -n);

% RCM reorder the jacobian
r = symrcm(B);
b = 2 * bandwidth(B(r,r)) + 1;

DIM.r = r;
DIM.b = b;

% Reorder Everything
DIM.PSISAT_AVG = PSISAT_AVG(r);
DIM.S_P = S_P(r);
DIM.S_E = S_E(r);
DIM.XZ = XZ(r, :);
DIM.NT = NT(r);
DIM.ST = ST(r, :);
DIM.DELTA = DELTA(r, :);
DIM.VOL = VOL(r, :);
end
