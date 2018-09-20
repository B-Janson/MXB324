close all
clear
clc
%% Part 0 Initialisation

T=linspace(0,24,25); %Units, maybe weeks?


[DIM] = GRIDCOORD;
[h, phi, k, S] = INITCOND(DIM);
%[M]=MATGEN(DIM,h);
%[S]=SGEN(DIM);


HEADVIS1(DIM, h)
HEADVIS2(DIM, phi)

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
