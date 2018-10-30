function BC = BOUNDARY_CONDITIONS(DIM, PARAMS, XZ, r_f, h)
%BOUNDARY_CONDITIONS function to return a vector of fluxes q . n
% Returned vector contains two elements and represents the same order as the
% DELTA matrix

x = XZ(1);
z = XZ(2);
BC = zeros(2, 1);

% Deal with LHS first
if x == 0
    if z == 0
        % Nothing because bedrock
    elseif z == DIM.HEIGHT
        % top left corner
        % top boundary is incoming rain
        BC(1) = -r_f;
        
        % left boundary should be river discharge but couldn't figure out
        % correct flux
        
%         if h + z > PARAMS.left_river(2)
%             BC(2) = PARAMS.K_r * 5 / 50;
%         end
    elseif z == PARAMS.left_river(1) % Bottom of river
        if h < 0 % Water table below river bed
            BC(2) = -PARAMS.K_r * (PARAMS.left_river(2) - PARAMS.left_river(1)) / PARAMS.left_river(3);
        else % water table above river bed
            BC(2) = -PARAMS.K_r * (PARAMS.left_river(2) - z - h) / PARAMS.left_river(3);
        end
    elseif PARAMS.left_river(1) < z && z < PARAMS.left_river(2)
        % in between river head and bed
        BC(1) = -PARAMS.K_r * (PARAMS.left_river(2) - (h + z)) / PARAMS.left_river(3);
        BC(2) = -PARAMS.K_r * (PARAMS.left_river(2) - (h + z)) / PARAMS.left_river(3);
    elseif z > PARAMS.left_river(2)
        % Above river head
        if h + z > PARAMS.left_river(2)
            BC(1) = -PARAMS.K_r * (PARAMS.left_river(2) - (h + z)) / PARAMS.left_river(3);
            BC(2) = -PARAMS.K_r * (PARAMS.left_river(2) - (h + z)) / PARAMS.left_river(3);
        end
    else
        % Nothing because bedrock
    end
elseif x == DIM.WIDTH % RHS
    if z == 0
        % Nothing because bedrock
    elseif z == DIM.HEIGHT
        % top right boundary
        BC(2) = -r_f;
        
        if h + z > PARAMS.right_river(2)
            BC(1) = PARAMS.K_r * 5 / 50;
        end
    elseif z == PARAMS.right_river(1) % Bottom of river
        if h < 0
            BC(2) = -PARAMS.K_r * 5 / 50;
        else
            BC(2) = -PARAMS.K_r * (5 - h) / 50;
        end
    elseif z > PARAMS.right_river(1)
        if h + z > PARAMS.right_river(2)
            BC(1) = PARAMS.K_r * 5 / 50;
            BC(1) = PARAMS.K_r * 5 / 50;
        else
            BC(1) = -PARAMS.K_r * h / 50;
            BC(2) = -PARAMS.K_r * h / 50;
        end
    else
        % Nothing because bedrock
    end
elseif z == 0 % Bottom
    % Nothing because bedrock
elseif z == DIM.HEIGHT
    BC = [-r_f; -r_f];
end

end
