function [x, m] = NEWTON_GMRES(A, b, x0, M, tol, max_krylov, diagnostics)

N = size(A, 1);
H = zeros(max_krylov + 1, max_krylov);
V = zeros(N, max_krylov + 1);

r = b - A * x0;
beta = norm(r, 2);
V(:, 1) = r / beta;
m = 0;
rnorm = inf;

while rnorm > beta * tol && m <= max_krylov
   m = m + 1;
   
   % Arnoldi (Modified Gram-Schmidt)
   V(:, m+1) = A * (M\V(:, m));
   for j = 1:m
       H(j, m) = V(:, j)' * V(:, m+1);
       V(:, m+1) = V(:, m+1) - H(j, m) * V(:, j);
   end
   H(m+1, m) = norm(V(:, m+1), 2);
   
   % Check for breakdown
   if abs(H(m+1, m)) < 1e-14
       fprintf('Invariant Krylov Subspace detected  at m=%d\n', m);
       y = H(1:m, 1:m) \ ([beta; zeros(m-1, 1)]);
       break;
   else
       V(:, m+1) = V(:, m+1) / H(m+1, m);
   end
   
   % Solve small m dimensional least squares problem for y
   rhs = [beta; zeros(m, 1)];
   y = H(1:m+1, 1:m) \ rhs;
   % Determine residual norm
   rnorm = norm(rhs - H(1:m+1, 1:m) * y);
   
   if diagnostics
       fprintf('m=%d, ||r_m||=%d tol=%d\n', m, rnorm, beta*tol);
   end
end

% Compute approximate solution
x = x0 + M\(V(:, 1:m) * y);


end