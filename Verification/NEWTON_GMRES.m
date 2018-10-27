function [x, m, rnorm] = NEWTON_GMRES(F,F_old,K, M, PARAMS, DIM, FVM_FUNC, h, h_old, S_old, phi_old, k_old, t)

max_krylov = PARAMS.gmres_max;
N = size(M, 1);
H = zeros(max_krylov + 1, max_krylov);
V = zeros(N, max_krylov + 1);

r = -F;
beta = norm(r, 2);
V(:, 1) = r / beta;
m = 0;

tol = PARAMS.gmres_tol;
x0 = zeros(N, 1);
[eta,Fk] = EW_ETA(K,F,F_old, PARAMS);
NORM=0;%get into loop
rnorm=inf;

%while  NORM < eta * Fk && m <= max_krylov   %Eisenstat Walker
while  rnorm > PARAMS.gmres_tol*beta && m <= max_krylov %Regular gmres
    m = m + 1;
    
    u = M \ V(:, m);
    
    if norm(h, 2) <= 1e-10
        delta = sqrt(eps);
    else
        delta = sqrt(eps) * norm(h, 2) / norm(u, 2);
    end
    
    JVm = (FVM_FUNC(DIM, h + delta*u, h_old, S_old, phi_old, k_old, t, PARAMS) - F) / delta;
    
    
    
    % Arnoldi (Modified Gram-Schmidt)
    V(:, m+1) = JVm;
    for j = 1:m
        H(j, m) = V(:, j)' * V(:, m+1);
        V(:, m+1) = V(:, m+1) - H(j, m) * V(:, j);
    end
    H(m+1, m) = norm(V(:, m+1), 2);
    rhs = [beta; zeros(m, 1)];
    
 %   NORM=norm(JVm-F,2);
    
    % Check for breakdown
    if abs(H(m+1, m)) < 1e-14
        fprintf('Invariant Krylov Subspace detected  at m=%d\n', m);
        Recomb=pinv(H(1:m+1,1:m));
        y=Recomb*rhs;
        break;
    else
        V(:, m+1) = V(:, m+1) / H(m+1, m);
    end
    
    %SVD method of calulating norm
    Recomb=pinv(H(1:m+1,1:m));
    y=Recomb*rhs;
    rnorm = norm(rhs - H(1:m+1, 1:m) * y);
    
    if ~PARAMS.debug
        fprintf('m=%d, ||r_m||=%d tol=%d\n', m, rnorm, beta*tol);
    end

    
end

    % Compute approximate solution
    x = x0 + M\(V(:, 1:m) * y);

PARAMS.eta_old=eta;


end