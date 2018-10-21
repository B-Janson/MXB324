function jacobian = JAC_FUNC(DIM, F, FVM_func, h, h_old, S_old, phi_old, k_old, t, PARAMS)
% Calculate the jacobian at the current time step

n = DIM.n;
m = DIM.m;

% Initialise the Matrix
jacobian = zeros(n*m, n*m);

% Is the solution the 0 vector?
a = norm(h, 2);

switch PARAMS.method
    case 'full'
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
            h_shift = h + step * shift_vec; % Create the shifted version of u for the finite difference

            % Evaluate the function with the shift
            Fhs = FVM_func(DIM, h_shift, h_old, S_old, phi_old, k_old, t, PARAMS);

            % Update the column of the jacobian with the finite difference formula
            jacobian(:,j) = (Fhs - F)/step;
        end
    case 'column'
        b = DIM.b;
        half_b = (b - 1) / 2;
        shifts = zeros(n*m, b); % Initialise array of shift vectors
        ep = zeros(1, b); % initialise vector of shift sizes
        Js = zeros(n*m, b);
        
        % Calculate shift vectors (vectors where bands do not overlap in the
        % Jacobian)
        for k = 1:n*m
           p = mod(k, b);
           if p == 0
               p = b;
           end
           shifts(k, p) = 1;
        end
        
        % calculate 2 norm of each shift vector
        for p = 1:b
            ep(p) = norm(shifts(:, p), 2);
        end

        % Calculate finite difference for each shift vector
        if a ~= 0
            ep = sqrt(eps) * a ./ ep; 
        else
            ep = sqrt(eps) ./ ep;
        end
        
        % Calculate the Jacobian at each shift vector
        for j = 1:b
            % shift h by vector
            h_shift = h + ep(j) * shifts(:,j); 
            % Evaluate the function with the shift
            Fhs = FVM_func(DIM, h_shift, h_old, S_old, phi_old, k_old, t, PARAMS);
            % Update the column of the jacobian with the finite difference formula
            Js(:,j) = (Fhs - F) / ep(j); % Store in array columns
        end

        % reorder shift vectors into the full jacobian
        % handle the first shift vectors up to the main diagonal separately
        for k = 1:b-1
           jacobian(1:k+half_b, k) = Js(1:k+half_b, k); 
        end
        
        % Do all remaining elements
        for k = 1:b
             idxs = find(shifts(b:n*m, k)) + b - 1; % find the ones in each shift vector after the ones previously accounted for
             for j = 1:length(idxs) % for each block left in the matrix
                 col_len_low = min(half_b, n*m - idxs(j)); % calculate if a block has the full length of the bandwidth or if it is shorter (near to the bottom of the Jacobian)
                 i = idxs(j)-half_b:idxs(j)+col_len_low; % calculate indexes of the block of elements
                 jacobian(i,idxs(j)) = Js(i,k); % store block in Js into corresponding column in the position Jacobian
             end
        end
end

jacobian = sparse(jacobian);

end
