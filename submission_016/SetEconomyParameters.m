function [] = SetEconomyParameters( beta_, eta_, alpha_, delta_, rho_, sigma_, N_, T_ )
   global beta eta alpha delta rho sigma N T previousC
   beta = beta_;
   eta = eta_;
   alpha = alpha_;
   delta = delta_;
   rho = rho_;
   sigma = sigma_;
   N = N_;
   T = T_;
   previousC = 0;
end
