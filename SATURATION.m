function [SAT] = SATURATION(DIM, h, i)



if h(i) < 0
    SAT = 0;
    VOL= DIM.VOL(i, :);
    for cell = 1:4
        %Pull from soil parameters table
        n = DIM.SP(DIM.ST(i,cell), 6);
        m = 1 - 1./n;
        alpha=DIM.SP(DIM.ST(i,cell), 5);
        SAT = SAT + VOL(cell) * (1 + (-alpha * h(i))^n)^-m;
    end
    %Average over control volume
    SAT = SAT / VOL(5);
else
    SAT = 1;
end

end