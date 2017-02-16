function [] = SetEconomyParameters( beta_, eta_, alpha_, delta_, rho_, sigma_, N_, T_ )

global beta eta alpha delta rho sigma N T est_ratio

beta = beta_;
eta = eta_;
alpha = alpha_;
delta = delta_;
rho = rho_;
sigma = sigma_;
N = N_;
T = T_;

K0 = 1;
Z0 = 1;

ratio_range = 0.1:0.005:1.0;
score = zeros(length(ratio_range), 1);
count = 1;
for ratio = ratio_range
    shocks = randn(T, N) * sigma;
    
    C = zeros(T+1, N);
    K = ones(T+1, N) * K0;
    Z = ones(T+1, N) * Z0;
    
    for t=1:T+1
        if t~=1
            Z(t,:) = exp(rho * log(Z(t-1,:)) + shocks(t-1,:));
        end
        W_t = (1 - delta) * K(t,:) + Z(t,:) .* (K(t,:) .^ alpha);
        if t~=T+1
            C(t,:) = ratio * W_t;
            K(t+1,:) = W_t - C(t,:);
        else
            C(t,:) = W_t;
        end
    end
    
    % Scoring
    eer_t = beta * ((C(1:end-1,:) ./ C(2:end,:)) .^ eta) .* (1 - delta + Z(2:end,:) * alpha .* (K(2:end,:) .^ (alpha - 1))) - 1;
    eer_t_avg = mean(eer_t, 2);
    eer = 1 ./ (1 + mean(abs(eer_t_avg)));
    score(count) = eer * 1000000;
    count = count + 1;
end
[~, idx] = max(score);
est_ratio = ratio_range(idx);
end