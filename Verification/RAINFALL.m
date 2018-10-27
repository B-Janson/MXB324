function [RF] = RAINFALL(PARAMS, t)
%RAINFALL Summary of this function goes here
%   Detailed explanation goes here

switch PARAMS.r_m
    case 1
        % constant rainfall and whether it's flood/drought/normal
        RF = PARAMS.r_f(PARAMS.r_t);
    case 2
        % cosine approximation of rainfall
        RF = PARAMS.r_f(PARAMS.r_t) * (cos(t * 2*pi / 365) + 1);
    case 3
        % realistic rainfall model
        %error('not yet implemented')
        RF = DALBY_RAIN(PARAMS.r_t, t);
end

end
