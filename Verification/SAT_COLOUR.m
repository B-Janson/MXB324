function sat_col = SAT_COLOUR()
%SAT_COLOUR Returns the colour scheme to plot the saturation contour
% Returns a 30 x 3 matrix of rgb values, where the r,g values decrease
% linearly.

sat_col = zeros(30, 3);

for i = 1:20
    sat_col(i, 1) = max(1.0 - 1.0 * (i-1) / 10, 0);
    sat_col(i, 2) = max(1.0 - 1.0 * (i-1) / 20, 0);
    sat_col(i, 3) = 1.0;
end

for i = 21:30
    sat_col(i, 3) = 1.0 - 0.5 * (i-21) / 10;
end

end
