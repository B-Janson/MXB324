function [DIM]=GRIDCOORD()
% GRIDCOORD returns the dimension of the grid

% Width and height of aquifer
W = 500;
H = 80;


DIM.x=linspace(0,W,3); %Keep it uniform for now
n=length(DIM.x);
DIM.n=n;
DIM.x=sort(DIM.x);


%The Discretisation we need in the y
DIM.z=linspace(0,H,3); %We Basic
m=length(DIM.z);
DIM.m=m;
DIM.z=sort(DIM.z);

NM=n*m;

%Create coordinate vector
[X,Z]=meshgrid(DIM.x,DIM.z);
X=X';
Z=Z';
XZ=[X(:),Z(:)];
DIM.XZ=XZ;


%% Create the distance matrix
%Create distance vectors dx=[L, R]
dx=zeros(DIM.n,2);
dx(1,2)=DIM.x(2)-DIM.x(1);
for i=2:DIM.n-1
    dx(i,1)=DIM.x(i)-DIM.x(i-1);
    dx(i,2)=DIM.x(i+1)-DIM.x(i);
end
dx(DIM.n,1)=DIM.x(DIM.n)-DIM.x(DIM.n-1);
%Create distance vectors dy=[D,U]
dz=zeros(DIM.m,2);
dz(1,2)=DIM.z(2)-DIM.z(1);
for i=2:DIM.m-1
    dz(i,1)=DIM.z(i)-DIM.z(i-1);
    dz(i,2)=DIM.z(i+1)-DIM.z(i);
end
dz(DIM.m,1)=DIM.z(DIM.m)-DIM.z(DIM.m-1);

%Create the big distances matrix
DELTA=zeros(m*n,4);

c=0;
for i=0:n:n*(m-1)
    DELTA(i+1:i+n,1:2)=dx;
    c=c+1;
    for j=1:n
        DELTA(i+j,3:4)=dz(c,:);
    end
end
%Each row [L,R,D,U]
DIM.DELTA=DELTA;

%Each row contains the volumes
%[UL,UR,DR,DL,TV]
DIM.VOL = zeros(n*m, 5);
for i = 1:m*n
    DIM.VOL(i, 1) = DELTA(i, 1) * DELTA(i, 4);%UL
    DIM.VOL(i, 2) = DELTA(i, 2) * DELTA(i, 4);%UR
    DIM.VOL(i, 3) = DELTA(i, 2) * DELTA(i, 3);%DR
    DIM.VOL(i, 4) = DELTA(i, 1) * DELTA(i, 3);%DL
    DIM.VOL(i,5)=sum(DIM.VOL(i,1:4));
end



% Set all cells to be Alluvium
%% Soil parameters (from PDF)
DIM.SP=[2.6,0.91,0.01,0.33,0.143,1.51;...
    0.08,0.0159,0.106,0.4686,1.04,1.3954;...
    3.9,1.17,0.0286,0.3658,2.8,2.239];
% Assign a node type to each vertex,
% Just the baby problem for now
% Add in the river nodes later

%Initialise the local soil/domain types
NT = zeros(NM, 1);
%Used to pull sub controll volume parameters from DIM.SP.
%1=alluvaim, 2=conf, 3=sandstone
%ordered [UL,UR,DR,DL]
ST = zeros(NM, 4); 


%Bottom L corner
NT(1) = 1;
ST(1,:)=[1,1,1,1];

%Bottom Row
NT(2:DIM.n-1) = 2;
for i=2:DIM.n-1
ST(i,:)=[1,1,1,1];
end

%Bottom R corner
NT(n)=3;
ST(n,:)=[1,1,1,1];

for i = n+1:NM-1
    if ((XZ(i,1) == 0) && (XZ(i,2) > 0) && (XZ(i,2) < 30))
        %Left Sandstone Boundary
        NT(i)=4;
    ST(i,:)=[1,1,1,1];
    elseif ( (XZ(i,1) > 0) && (XZ(i,1) < 500) && (XZ(i,2) > 0) && (XZ(i,2) < 30) ) || ...
            ( (XZ(i,1) > 350) && (XZ(i,1) < 500) && (XZ(i,2) > 0) && (XZ(i,2) < 80) )
        
        % Sandstone Interior
        % Pump at (450, 10)
        
        if ((XZ(i,1) == 450) && (XZ(i,2) == 10))
            NT(i) = 5;
            % Evapotranspiration zone (z >= 76)
        ST(i,:)=[1,1,1,1];
        elseif XZ(i,2) >= 76
            NT(i)=6;
            ST(i,:)=[1,1,1,1];
            % Normal sandstone interior point
        else
            NT(i)=7;
        ST(i,:)=[1,1,1,1];
        end
    elseif ((XZ(i,1) == 500) && (XZ(i,2) > 0) && (XZ(i,2) < 80))
        % Right sandstone boundary
        if XZ(i,2) < 76
            NT(i)=8;
            ST(i,:)=[1,1,1,1];
            % Evapotranspiration zone
        else
            NT(i)=9;
            ST(i,:)=[1,1,1,1];
        end
    elseif (XZ(i,2) == 30)
        % Sandstone/Alluvium & Boundary
        if (XZ(i,1) == 0)
            NT(i)=10;
            ST(i,:)=[1,1,1,1];
        elseif ((XZ(i,1) > 0) && (XZ(i,1) < 50))
            % Sandstone/Alluvium Interface
            NT(i)=11;
            ST(i,:)=[1,1,1,1];
        elseif (XZ(i,1) == 50)
            % Sandstone/Alluvium/Confining interface
            NT(i)=12;
            ST(i,:)=[1,1,1,1];
        elseif ((XZ(i,1) > 50) && (XZ(i,1) < 350))
            % Sandstone/Confining horizontal interface
            NT(i)=13;
            ST(i,:)=[1,1,1,1];
        elseif (XZ(i,1) == 350)
            % Sandstone/Confining corner
            NT(i)=14;
            ST(i,:)=[1,1,1,1];
        end
    elseif ((XZ(i,1) == 0) && (XZ(i,2) > 30) && (XZ(i,2) < 80))
        % Alluvium left boundary
        if XZ(i,2) < 78
            NT(i)=15;
            ST(i,:)=[1,1,1,1];
        else
            % Alluvium boundary evapotranspiration
            NT(i)=16;
            ST(i,:)=[1,1,1,1];
        end
    elseif ( (XZ(i,1) > 0) && (XZ(i,1) < 50) && (XZ(i,2) > 30) && (XZ(i,2) < 80) ) || ...
            ( (XZ(i,1) > 0) && (XZ(i,1) < 350) && (XZ(i,2) > 40) && (XZ(i,2) < 80) )
        % Alluvium Interior
        if ((XZ(i,1) == 100) && (XZ(i,2) == 50))
            % Pump at (100, 50)
            NT(i)=17;
            ST(i,:)=[1,1,1,1];
        elseif XZ(i,2) >= 78
            % Evapotranspiration
            NT(i)=18;
            ST(i,:)=[1,1,1,1];
        else
            NT(i)=19;
            ST(i,:)=[1,1,1,1];
        end
    elseif ((XZ(i,1) == 50) && (XZ(i,2) > 30) && (XZ(i,2) < 40))
        % Alluvium/Confining interface
        NT(i)=20;
        ST(i,:)=[1,1,1,1];
    elseif ((XZ(i,1) > 50) && (XZ(i,1) < 350) && (XZ(i,2) > 30) && (XZ(i,2) < 40))
        % Confining interior
        NT(i)=21;
        ST(i,:)=[1,1,1,1];
    elseif ((XZ(i,1) == 350) && (XZ(i,2) > 30) && (XZ(i,2) < 40))
        % Confining/Sandstone vertical interface
        NT(i)=22;
        ST(i,:)=[1,1,1,1];
    elseif ((XZ(i,1) == 50) && (XZ(i,2) == 40))
        % Confining/Alluvium corner
        NT(i)=23;
        ST(i,:)=[1,1,1,1];
    elseif ((XZ(i,1) > 50) && (XZ(i,1) < 350) && (XZ(i,2) == 40))
        % Confining/Alluvium horizontal interface
        NT(i)=24;
        ST(i,:)=[1,1,1,1];
    elseif ((XZ(i,1) == 350) && (XZ(i,2) == 40))
        % Confining/Alluvium/Sandstone upper right corner
        NT(i)=25;
        ST(i,:)=[1,1,1,1];
    elseif ((XZ(i,1) == 350) && (XZ(i,2) > 40) && (XZ(i,2) < 80))
        % Alluvium/Sandstone interface
        if XZ(i,2) < 76
            NT(i)=26;
            ST(i,:)=[1,1,1,1];
        elseif XZ(i,2) < 78
            % First half of evapo
            NT(i)=27;
            ST(i,:)=[1,1,1,1];
        else
            % Second half of evapo
            NT(i)=28;
            ST(i,:)=[1,1,1,1];
        end
    elseif 1 == 0 %((XY(i,1) == 0) && ((XY(i,2) == 60))
        %Bed/River
        NT(i)=28;
        ST(i,:)=[1,1,1,1];
    elseif 1 == 0 %((XY(i,1) == 0) && ((XY(i,2) > 60))) && ((XY(i,2) <80))
        %River Boundary
        
        
    elseif ((XZ(i,1) == 0) && (XZ(i,2) == 80))
        % Top left node point River/Rain, L bedrock rain for baby problem
        NT(i)=29;
        ST(i,:)=[1,1,1,1];
    elseif ((XZ(i,1) > 0) && (XZ(i,1) < 350) && (XZ(i,2) == 80))
        % Alluvium/Rain boundary
        NT(i)=30;
        ST(i,:)=[1,1,1,1];
    elseif ((XZ(i,1) == 350) && (XZ(i,2) == 80))
        % Alluvium/Sandstone/Rain Boundary
        NT(i)=31;
        ST(i,:)=[1,1,1,1];
    elseif ((XZ(i,1) > 350) && (XZ(i,1) < 500) && (XZ(i,2) == 80))
        % Sandstone/Rain Boundary
        NT(i)=32;
        ST(i,:)=[1,1,1,1];
    end
end

NT(NM)=33;
ST(NM,:)=[1,1,1,1];

B=gallery('tridiag',NM,1,1,1);
L=n;
U=NM;
for i = 1:n*(m-1)
    B(L,i)=1;
    B(U-n,U)=1;
    L=L+1;
    U=U-1;
end

figure()
spy(B)


r=symrcm(B);
Weightloss=2*(bandwidth(B)-bandwidth(B(r,r)))

figure()
spy(B(r,r))

DIM.r=r;
DIM.b=2*b2+1

DIM.XZ=DIM.XZ(r,:);
DIM.NT=NT(r);
DIM.ST=ST(r,:);
end