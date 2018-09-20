function [DIM]=GRIDCOORD()


%The Discretisation we need in the y coordinate

DIM.x=linspace(0,500,51); %Keep it uniform for now
n=length(DIM.x);
DIM.n=n;
DIM.x=sort(DIM.x);
DIM.dx=zeros(1,DIM.n-1);
for i=1:DIM.n-1
    DIM.dx(i)=DIM.x(i+1)-DIM.dx(i);
end

%The Discretisation we need in the y
DIM.y=linspace(0,80,6); %We Basic
m=length(DIM.y);
DIM.m=m;
DIM.y=sort(DIM.y);
DIM.dy=zeros(1,DIM.m-1);


%Create coordinate vector
[X,Y]=meshgrid(DIM.x,DIM.y);
X=X';
Y=Y';
XY=[X(:),Y(:)];
DIM.XY=XY;

%Assign a node type to each vertex,
%Just the baby problem for now
%Add in the river nodes later
NT=zeros(n*m,1);

%Bottom L corner
NT(1)=1;
%Bottom Row
NT(2:DIM.n-1)=2;
%Bottom R corner
NT(n)=3;
for i=n+1:n*m-1
    if ((XY(i,1) == 0) && ((XY(i,2) > 0))) && ((XY(i,2) <80))
        %Left Sandstone Boundary
        NT(i)=4;
    elseif (((XY(i,1) > 0) && (XY(i,1)) < 500) && ((XY(i,2) > 0)) && ((XY(i,2) <30))) || (((XY(i,1) > 350) && (XY(i,1)) < 500) && (XY(i,2) > 0) && (XY(i,2) <80))
        %Sandstone Interior
        if (~(XY(i,1) == 450) && (XY(i,2) == 10))
            NT(i)=5;
        elseif XY(i,2) >= 76
            NT(i)=6;
        else
            NT(i)=7;
        end
    elseif ((XY(i,1) == 500) && ((XY(i,2) > 0))) && ((XY(i,2) <80))
        %Right sandstone boundary
        if XY(i,2) < 76
            NT(i)=8;
        else
            NT(i)=9;
        end
    elseif (XY(i,2) == 30) && (XY(i,1) == 0)
        %S/A Boundary
        NT(i)=10;
    elseif ((XY(i,2) == 30) && ((XY(i,1) > 0))) && ((XY(i,1) < 50))
        %S/A Interface
        NT(i)=11;
    elseif (XY(i,2) == 30) && (XY(i,1) == 50)
        %S/A/C lower corner
        NT(i)=12;
    elseif ((XY(i,2) == 30) && ((XY(i,2) > 0))) && ((XY(i,2) <80))
        %S/C horizontal interface
        NT(i)=13;
    elseif (XY(i,2) == 30) && (XY(i,1) == 350)
        %S/C corner
        NT(i)=14;
    elseif (XY(i,1) == 0) && (XY(i,2) > 30)
        %A left boundary
        if XY(i,2) < 78
            NT(i)=15;
        else
            NT(i)=16;
        end
    elseif (((XY(i,1) > 0) && (XY(i,1)) < 50) && ((XY(i,2) > 30)) && ((XY(i,2) > 80))) || (((XY(i,1) > 0) && (XY(i,1)) < 500) && (XY(i,2) > 40) && (XY(i,2) <80))
        %A Interior
        if (~(XY(i,1) == 100) && (XY(i,2) == 50))
            NT(i)=17;
        elseif XY(i,2) >=78
            NT(i)=18;
        else
            NT(i)=19;
        end
    elseif ((XY(i,1) == 50) && ((XY(i,2) > 30))) && ((XY(i,2) < 40))
        %A/C interface
        NT(i)=20;
    elseif (((XY(i,1) > 50) && (XY(i,1)) < 350) && ((XY(i,2) > 30))) && ((XY(i,2) < 40))
        %C interior
        NT(i)=21;
    elseif ((XY(i,1) == 350) && ((XY(i,2) > 30))) && ((XY(i,2) <50))
        %C/S vertical interface
        NT(i)=22;
    elseif ((XY(i,1) == 50)) && ((XY(i,2) == 40))
        %C/A corner
        NT(i)=23;
    elseif ((XY(i,2) == 40) && ((XY(i,1) > 50))) && ((XY(i,2) < 350))
        %C/A horizontal interface
        NT(i)=24;
    elseif ((XY(i,1) == 350)) && ((XY(i,2) == 40))
        %C/A/S upper corner
        NT(i)=25;
    elseif ((XY(i,1) == 350) && ((XY(i,2) > 40))) && ((XY(i,2) < 80))
        %A/S interface
        if XY(i,2) < 78
            NT(i)=26;
        else
            NT(i)=27;
        end
    elseif 1 == 0 %((XY(i,1) == 0) && ((XY(i,2) == 60))
        %Bed/River
        NT(i)=28;
    elseif 1 == 0 %((XY(i,1) == 0) && ((XY(i,2) > 60))) && ((XY(i,2) <80))
        %River Boundary
        
        
    elseif ((XY(i,1) == 0)) && ((XY(i,2) == 80))
        %River/Rain, L bedrock rain for baby problem
        NT(i)=29;
    elseif ((XY(i,2) == 80) && ((XY(i,1) > 0))) && ((XY(i,1) < 350))
        %A/Rain boundary
        NT(i)=30;
    elseif ((XY(i,1) == 350)) && ((XY(i,2) == 80))
        %A/S/Rain Boundary
        NT(i)=31;
    elseif ((XY(i,2) == 80) && ((XY(i,1) > 350))) && ((XY(i,1) <500))
        %S/Rain Boundary
        NT(i)=32;
    elseif (XY(i,1) == 500) && (XY(i,2) == 80)
        %S/bedrock/rain
        NT(i)=33;
    end
end
NT(DIM.n*DIM.m)=34;

B=gallery('tridiag',n*m,1,1,1);
L=n;
U=n*m;
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