function jacobian = JAC_FUNC(DIM, F, FVM_func, h, h_old, S_old, phi_old, k_old, t, PARAMS)

n = DIM.n;
m = DIM.m;

% Initialise the Matrix
jacobian = zeros(n*m, n*m);

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
    h_shift = h + step*shift_vec; % Create the shifted version of u for the finite difference
    % Evaluate the function with the shift
    Fhs = FVM_func(DIM, h_shift, h_old, S_old, phi_old, k_old, t, PARAMS);
    % Update the column of the jacobian with the finite difference formula
    jacobian(:,j) = (Fhs - F)/step;
end

end
