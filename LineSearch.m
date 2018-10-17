function h_star = LineSearch(DIM, FVM_func, h, del_h, h_old, S_old, phi_old, k_old, t, PARAMS)

lambda = 1;
min_lambda = 10^-6;
h_star = h + lambda * del_h;

base_error = norm(FVM_func(DIM, h, h_old, S_old, phi_old, k_old, t, PARAMS));
current_error = norm(FVM_TEST(DIM, h_star, h_old, S_old, phi_old, k_old, t, PARAMS));

% fprintf('Line searching ENABLED')

while current_error > base_error && lambda > min_lambda
    lambda = lambda/2;
    h_star = h + lambda * del_h;
    current_error = norm(FVM_func(DIM, h_star, h_old, S_old, phi_old, k_old, t, PARAMS));
%     fprintf("Line search ENABLED \n")
end

if lambda < min_lambda
    fprintf("min lambda reached\n")
end
