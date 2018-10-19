function phi_true = PHI_TRUE(DIM, PARAMS, t)
%PHI_TRUE Calculate the analytic solution for water content

phi_true = PARAMS.r_f(1) * t / DIM.HEIGHT + 0.0767;
end
