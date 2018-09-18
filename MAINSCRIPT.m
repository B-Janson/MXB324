close all
clear
%% Part 0 Initialisation

T=linspace(0,24,25)%Units, maybe weeks?


[DIM]=GRIDCOORD;
U0=INITCOND(DIM);
%[M]=MATGEN(DIM,h);
%[S]=SGEN(DIM);


HEADVIS1(DIM,U0)
PHI=WCONT(DIM,U0);
HEADVIS2(DIM,PHI)

%% PART 1 Convergence?







%Storage
%U is an (m*n)x1 matrix stroing the initial condition
%
%
%for start to finish
%u(:,end+1)=last solution
%Solve system
%T(i)=timesteptaken this step
%end
