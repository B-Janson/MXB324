function J = jac_funcOG(DIM, F, FVM_func, h, h_old, S_old, phi_old, k_old, t, PARAMS)

n = DIM.n;
m = DIM.m;


if strcmp(PARAMS.method,'full')
    % Initialise the Matrix
    J = zeros(n*m, n*m);
    
    a = norm(h, 2); % Is the solution the 0 vector?
    
    % Find a suitable finite difference
    if a ~= 0
        step = sqrt(eps)*a;
    else
        step = sqrt(eps);
    end
    % eps - machine epsilon or the error of the floating point arithmetic in
    % MATLAB
    for j = 1:n*m
        % create the shift vector and put a 1 corresponding to the column of
        % the Jacobian being calculated
        shift_vec = zeros(n*m,1);
        shift_vec(j) = 1;
        nudge = h + step*shift_vec; % Create the shifted version of u for the finite difference
        %0 Evaluate the function with the shift
        Fnudge = FVM_func(DIM, nudge, h_old, S_old, phi_old, k_old, t, PARAMS);
        % Update the column of the jacobian with the finite difference formula
        J(:,j) = (Fnudge - F)/step;
    end
    [~,c]=chol(J);
    if c ~= 0
        error('Not SPD');
    end
elseif strcmp(PARAMS.method,'column')
    %Creates the column based finite differences matrix
    b=DIM.b;
    B=2*b+1;
    S=zeros(n*m,b);
    Jcol=S;
    %Finds the s-vectors
    for i=1:b:n*m
        for j=1:b
            S(i+j-1,j)=1;
        end
    end
    %Generate jacobian
    step=S;
    
    for i=1:b
        a=norm(h,2);
        s=norm(S(:,i),2);
        % Find a suitable finite difference
        if a ~= 0
            step(:,i) = step(:,i)*sqrt(eps)*a/s;
        else
            step(:,i) = step(:,i)*sqrt(eps)/s;
        end
        Fnudge=FVM_func(DIM, h+step(i,:), h_old, S_old, phi_old, k_old, t, PARAMS);
        Jcol(:,i)=(Fnudge-F)./step(:,i);
    end
    
    %Put elements into jacobian
    
    J=zeros(m*n,n*m);
    
    %First Bands
    
    for i=1:b-2
        J(i,1:1+i)=Jcol(i,1:i+1);      
    end    
    %Middle Bands
    
    start=0;    
    for i=b-1:m*n-1
        J(i,1+start:b+start)=Jcol(i,:);
        start=start+1;        
    end
    %Last Band
    J(m*n,m*n-b+2:m*n)=Jcol(m*n,2:b);
    
    [~,c]=chol(J);
    if c ~= 0
        eig(J)
        error('Not SPD');
    end
end
end


