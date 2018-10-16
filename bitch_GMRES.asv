function [J]=bitch_GMRES(J,J_old,F_old,DIM,PARAMS)

tol=PARAMS.gmres_tol; %max residual
max=PARAMS.gmres_max; %max iterations

n=DIM.n;
m=DIM.m;

x0=zeroes(m*n,1);

Beta=norm(b,2);

V1=b/Beta;


C=1;

H=zeros(C,C);
V=V1;

Beta_e1=[Beta;zeros(n*m-1,1)];

while (r > tol) || (C <= max)


    %Start Arnoldi Algorithm
    V(:,C+1)=A*V(:,m);
    for i=1:C
        H(i,C)=V(:,i)'*V(:,C+1);
        V(:,C+1)=V(:,m+1)-H(i,C)*V(:,i);
    end
    
    H(C+1,C)=norm(V(:,C+1),2);
    if H(m+1,m) ~= 0
       V(:,C+1)=V(:,C+1)/H(C+1,C); 
    end       
    
    %Start Least Squares


    C=C+1;
end

end