function h_star = LineSearch(DIM, FVM_func, h, del_h, h_old, S_old, phi_old, k_old, t, PARAMS)
% LINESEARCH  returns a line searched value for h
%   h_star = h + lambda * del_h, where lambda gets halved until a better
%   solution is found

% initialise lambda and minimum value
lambda = 1;
min_lambda = 10^-3;
% initialise alpha constant
alpha = 1e-4;
% initialise sigma constant
sigma = 0.5;

% calculate base error
h_star = h + lambda * del_h;

g_base = norm(FVM_func(DIM, h, h_old, S_old, phi_old, k_old,  t, PARAMS), 2)^2;

g_star = norm(FVM_func(DIM, h_star, h_old, S_old, phi_old, k_old, t, PARAMS), 2)^2;

% while we haven't reduced error sufficiently
while g_star >= (1 - 2 * alpha * lambda) * g_base && lambda > min_lambda
    % Halve lambda
    lambda = sigma * lambda;
    % Get a new value for h
    h_star = h + lambda * del_h;
    % Update this error
    g_star = norm(FVM_func(DIM, h_star, h_old, S_old, phi_old, k_old, t, PARAMS), 2)^2;
end

end
