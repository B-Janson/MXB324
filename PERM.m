function [k] = PERM(DIM, h, S, i)




if h(i) < 0
    k = 0;
    %Find volumes
    VOL = DIM.VOL(i, :);
    for cell = 1:4
        %Pull parameters from table
        n = DIM.SP(DIM.ST(i,cell), 6);
        m = 1 - 1./n;
        %Calculate the permiability in each control volume
        k = k + VOL(cell) * (sqrt(S(i)) * (1 - (1 - S(i)^(1/m))^m)^2);
    end
    %Average over the cell
    k = k /VOL(5);
else
    k = 1;
end

end