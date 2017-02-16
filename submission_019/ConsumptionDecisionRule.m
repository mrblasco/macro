function [ C_t ] = ConsumptionDecisionRule( K_t, Z_t )

   global beta eta alpha delta rho sigma N T est_ratio

   W_t = (1 - delta) * K_t + Z_t .* (K_t .^ alpha);
   C_t = est_ratio * W_t;

end