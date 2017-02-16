function [ C_t ] = ConsumptionDecisionRule(K_t, Z_t)
   global beta eta alpha delta rho sigma N T previousC
   W_t = (1 - delta) * K_t + Z_t .* (K_t .^ alpha);
   
   %Baseline
   a=1/beta;
   c=times(Z_t, alpha);
   d=power(K_t, alpha-1);
   e=times(c,d);
   f=1-delta+e;
   g=ldivide(a,f);
   h=power(g, (1.0/eta));
   j=ldivide(1,h);
   
   if(previousC==0)
       C_t=times(0.0001, W_t);
   else
       C_t=times(j, previousC);
   end
   
   for i=1:length(C_t)
       if(C_t(i)>W_t(i))
          C_t(i)=times(0.99, W_t(i));
       elseif(C_t(i)<=0)
          C_t(i)=times(0.01, W_t(i));
       end
   end
      
   previousC=C_t;
end

