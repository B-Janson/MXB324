function [PHI]=WATER_CONTENT(DIM, h, S, i)
% Water content function

PHI = 0;
%Pull volumes of the cell
VOL = DIM.VOL(i, :);

if h(i) < 0 %Unsaturated
    for cell = 1:4
        %Pull parameters from table
        PR =DIM.SP(DIM.ST(i,cell), 3);
        PS =DIM.SP(DIM.ST(i,cell), 4);
        %Find water content of the sub control volume
        PHI = PHI + VOL(cell) * (PR + S(i) * (PS - PR));
    end
    %Average over whole cell
    PHI = PHI / VOL(5);
else %Saturated
    for cell = 1:4
        %Pull parameter from table
        PS =DIM.SP(DIM.ST(i,cell), 4);
        %Find water content of the sub control volume
        PHI = PHI + VOL(cell) * PS;
    end
    %Average over whole cell
    PHI = PHI / VOL(5);
end

end