function BC = BOUNDARY_CONDITIONS(PARAMS)

vert_lines = (PARAMS.m - 1) * 2;
horiz_lines = (PARAMS.n - 1) * 2;

% River locations (either left or right of domain).
% If no river wanted, set river_left or river_right = [0, 0];

% Z location of left river [head, bottom]
river_left = [65, 60];

if river_left(1) < river_left(2)
   error('the river head must be above the river bottom\')
end

% Z location of right river [head, bottom]
river_right = [0, 0];

if river_right(1) < river_right(2)
   error('the river head must be above the river bottom\')
end

% General left boundary condition
left_BC = zeros(1, vert_lines);

% General right boundary condition
right_BC = zeros(1, vert_lines);

% General bottom boundary conditions
bottom_BC = zeros(1, horiz_lines);

BC.river_left = river_left;
BC.river_right = river_right;

BC.left = left_BC;
BC.right = right_BC;
BC.bottom = bottom_BC;


end