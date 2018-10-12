function J = jac_func(DIM, F, FVM_func, h, h_old, S_old, phi_old, k_old, t, PARAMS)

n = DIM.n;
m = DIM.m;


if strcmp(PARAMS.method,'Full')
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
elseif strcmp(PARAMS.method,'Column')
    %Creates the column based finite differences matrix
    b=DIM.b;
    B=2*b+1;
    S=zeros(n*m,B);
    Jcol=S;
    step=ones(n*m,B);
    %Finds the s-vectors
    for i=1:B:n*m
        for j=1:B
            S(i+j-1,j)=1;
        end
    end
    %Generate jacobian
    
    for i=1:B
        
        a=norm(h,2);
        s=norm(S(:,i),2);
        % Find a suitable finite difference
        if a ~= 0
            step(:,i) = step(:,i)*sqrt(eps)*a/s;
        else
            step(:,i) = step(:,i)*sqrt(eps)/s;
        end
        Fnudge=FVM_func(DIM, h+step(i,:), h_old, S_old, phi_old, k_old, t, PARAMS);
        Jcol(:,i)=(Fnudge-F)./step(:,i)   ;
    end
    
    %Put elements into jacobian
    
    J=zeros(m*n,n*m);
    
    %First Band
    fin=b;
    for i=1:b
        J(i,1:1+fin)=Jcol(i,1:1+fin);
        fin=fin+1;
        
    end
    %Middle Bands
    for i=b+1:m*n-b-1
        J(i,i-b:i+b)=Jcol(i,:);
        
    end
    %Last Band
    begin=1;
    for i=m*n-b:m*n
        J(i,begin+b-1:m*n)=Jcol(i,begin:B);
        begin=begin+1;
    end   
    [~,c]=chol(J);
    if c ~= 0
        error('Not SPD');
        
        eig(J)
    end
end
end


