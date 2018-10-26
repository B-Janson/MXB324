function [eta] = EW_ETA(k,eta_old,F,F_old, PARAMS)
%Need to precalculate ||F(0)||_inf


if (k == 0)
    eta=PARAMS.eta_max;
else
    Fk=norm(F,inf);
    Fk_old=nrom(F_old,inf);
    eta_r=PARAMS.lambda*(Fk/Fk_old)^PARAMS.alpha;
    C=PARAMS.lambda*(eta_old)^PARAMS.alpha;
    if  C <= 0.1
        eta_s=min(PARAMS.eta_max,eta_r);
    else
        eta_s=min(PARAMS.eta_max,max(eta_r,C));
    end
    eta=min(PARAMS.eta_max,max(eta_s,0.5*(PARAMS.TA+PARAMS.TR*PARAMS.F0)/Fk));
end

end
