function [ C_t ] = ConsumptionDecisionRule( K_t, Z_t )

   global beta eta alpha delta rho sigma N T

   W_t = (1 - delta) * K_t + Z_t .* (K_t .^ alpha);
   
   %linear LQ Approximation
   Kstar = (beta*alpha/(1-beta))^(1/(1-alpha));

   Cstar = (1-beta) * Kstar / (beta*alpha);
   phi = 1 + 1/beta + (1-alpha) * (1-beta) * Cstar/(Kstar*eta);
   q = beta * ((1-rho)*Cstar+rho*(1-beta)*Cstar/eta);

   lamda = phi/2 - sqrt((phi*beta)^2 - 4*beta)/(2*beta);
   if abs(lamda)<=1
      K_t_1 = (1-lamda)*Kstar + lamda*K_t + q*lamda*log(Z_t)/(1-beta*rho*lamda);
      C_t = W_t - K_t_1;
   end
   
   %Exception handle
   if any(C_t>W_t) || any(C_t<=0)
        C_t = 0.5 * W_t;
   end
end

