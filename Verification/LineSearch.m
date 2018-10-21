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

% the error using the current h
% base_error = norm(FVM_func(DIM, h, h_old, S_old, phi_old, k_old, t, PARAMS));

% the error using the h_star
% current_error = norm(FVM_func(DIM, h_star, h_old, S_old, phi_old, k_old, t, PARAMS));

g_base = norm(FVM_func(DIM, h, h_old, S_old, phi_old, k_old,  t, PARAMS), 2)^2;

g_star = norm(FVM_func(DIM, h_star, h_old, S_old, phi_old, k_old, t, PARAMS), 2)^2;

% while we haven't reduced error
while g_star >= (1 - 2 * alpha * lambda) * g_base && lambda > min_lambda
    lambda = sigma * lambda;
    h_star = h + lambda * del_h;
    g_star = norm(FVM_func(DIM, h_star, h_old, S_old, phi_old, k_old, t, PARAMS), 2)^2;
end

if lambda < min_lambda && PARAMS.debug == true
    fprintf("min lambda reached\n")
end

end
