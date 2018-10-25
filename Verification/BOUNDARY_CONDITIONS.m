function BC = BOUNDARY_CONDITIONS(DIM, PARAMS, XZ, t, h)

x = XZ(1);
z = XZ(2);
BC = zeros(2, 1);

switch PARAMS.r_m
    case 1
        r_f = PARAMS.r_f(PARAMS.r_t);
    case 2
        r_f = PARAMS.r_f(PARAMS.r_t) * (cos(t * 2*pi / 365) + 1);
    case 3
        error('not yet implemented');
end

% Deal with LHS first
if x == 0
    if z == 0
        % Nothing because bedrock
    elseif z == DIM.HEIGHT
        BC(1) = -r_f;
        
%         if h + z > PARAMS.left_river(2)
%             BC(2) = PARAMS.K_r * 5 / 50;
%         end
    elseif z == PARAMS.left_river(1) % Bottom of river
        if h < 0
            BC(2) = -PARAMS.K_r * (PARAMS.left_river(2) - PARAMS.left_river(1)) / PARAMS.left_river(3);
        else
            BC(2) = -PARAMS.K_r * (PARAMS.left_river(2) - z - h) / PARAMS.left_river(3);
        end
    elseif PARAMS.left_river(1) < z && z < PARAMS.left_river(2)
        BC(1) = -PARAMS.K_r * (PARAMS.left_river(2) - (h + z)) / PARAMS.left_river(3);
        BC(2) = -PARAMS.K_r * (PARAMS.left_river(2) - (h + z)) / PARAMS.left_river(3);
%         if h + z < PARAMS.left_river(1)
%             BC(1) = -PARAMS.K_r * (PARAMS.left_river(2) - PARAMS.left_river(1)) / PARAMS.left_river(3);
%             BC(2) = -PARAMS.K_r * (PARAMS.left_river(2) - PARAMS.left_river(1)) / PARAMS.left_river(3);
%         elseif h + z <= PARAMS.left_river(2)
%             BC(1) = PARAMS.K_r * (PARAMS.left_river(2) - (h + z)) / PARAMS.left_river(3);
%             BC(2) = PARAMS.K_r * (PARAMS.left_river(2) - (h + z)) / PARAMS.left_river(3);
%         end
    elseif z > PARAMS.left_river(2)
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
