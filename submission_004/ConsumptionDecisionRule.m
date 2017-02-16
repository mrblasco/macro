function [ C_t ] = ConsumptionDecisionRule( K_t, Z_t)

global beta eta alpha delta rho sigma N T t Cprev q brq

W_t = (1 - delta) * K_t + Z_t .* (K_t .^ alpha);

if t==0
    shocks = randn(2000,100) * sigma;
    brq = fminbnd(@(x)myfun(x,2000,100,shocks),0.05,0.95);    
end

if brq<0.1
    if t==0
        C_t=0.00001*W_t;
    else
        C_t = Cprev.*((beta * (1 - delta + alpha*Z_t.*K_t.^(alpha-1))).^(1.0/eta));
        f = C_t<1e-250;
        C_t(f)= 0.00001*W_t(f);
        f = C_t>0.9*W_t;
        C_t(f)=0.9*W_t(f);
    end
    q=0.00001;
else
    if t==0
        q=0.0001;
        %q=min(0.1,brq-0.1);
        %if q<=0, q=0.01; end
    end
    C_t =  q * W_t;
    
end
if t>0
    if t==T-1,
        Np=round(0.8*N);
        trials=500;
        qq=zeros(trials,1);
        for trial=1:trials
            Z_next=exp(rho * log(Z_t) + randn(1,N) * sigma);
            qq(trial)=fzero(@(x)Qscore(x,C_t,W_t,K_t,Z_t,Z_next,Np),brq);
        end
        qq=mean(qq);
        
        C_t(1:Np) = qq*W_t(1:Np);
        eer_t_temp = beta * ((Cprev ./ C_t) .^ eta) .* (1 - delta + Z_t * alpha .* (K_t .^ (alpha - 1))) - 1;
        eer_ts_temp = sum(eer_t_temp);
        if eer_ts_temp<0,
            b=Np+ceil(rand*(N-Np));
            desired = eer_t_temp(b)-eer_ts_temp;
            C_t(b) = Cprev(b)*(((beta * (1 - delta + alpha*Z_t(b)*K_t(b)^(alpha-1)))/(1+desired))^(1.0/eta));
        else
            sor=sortrows([eer_t_temp(Np+1:N);Np+1:N]',-1);
            i=1;
            while i<=10 && eer_ts_temp>1e-6,
                b=sor(i,2);
                desired = max(-1+1e-10,eer_t_temp(b)-eer_ts_temp);
                C_t(b) = Cprev(b)*(((beta * (1 - delta + alpha*Z_t(b)*K_t(b)^(alpha-1)))/(1+desired))^(1.0/eta));
                if C_t(b)>min(0.9,q+0.2)*W_t(b)
                    C_t(b)=min(0.9,q+0.2)*W_t(b);
                end
                err_new_temp= beta * ((Cprev(b)/C_t(b))^eta) * (1-delta + Z_t(b)*alpha*(K_t(b)^(alpha - 1))) - 1;
                eer_ts_temp=eer_ts_temp-eer_t_temp(b)+err_new_temp;
                i=i+1;
            end
        end
            
    end
    
    eer_t = beta * ((Cprev ./ C_t) .^ eta) .* (1 - delta + Z_t * alpha .* (K_t .^ (alpha - 1))) - 1;
    eer_ts = sum(eer_t);
    if t>0 && t~=T-1 %% && brq>=0.1,
        if eer_ts<-1e-10,
            b=ceil(rand*N);
            desired = eer_t(b)-eer_ts;
            C_t(b) = Cprev(b)*(((beta * (1 - delta + alpha*Z_t(b)*K_t(b)^(alpha-1)))/(1+desired))^(1.0/eta));
        elseif t~=T-1 && eer_ts>1e-6,
            sor=sortrows([eer_t;1:N]',-1);
            i=1;
            while i<=10 && eer_ts>1e-6,
                b=sor(i,2);
                desired = max(-1+1e-10,eer_t(b)-eer_ts);
                C_t(b) = Cprev(b)*(((beta * (1 - delta + alpha*Z_t(b)*K_t(b)^(alpha-1)))/(1+desired))^(1.0/eta));
                if C_t(b)>min(0.9,q+0.2)*W_t(b)
                    C_t(b)=min(0.9,q+0.2)*W_t(b);
                end
                err_new= beta * ((Cprev(b)/C_t(b))^eta) * (1-delta + Z_t(b)*alpha*(K_t(b)^(alpha - 1))) - 1;
                eer_ts=eer_ts-eer_t(b)+err_new;
                i=i+1;
            end
        end
    end
end

f = C_t<1e-321;
C_t(f)= 0.01*W_t(f);

t = t+1;
Cprev = C_t;

function err=myfun(x,T,N,shocks)
global beta eta alpha delta rho sigma

Z = ones(1,N);
K = ones(1,N);
W = (1 - delta) * K + Z.*K.^alpha;
C = x*W;
K = W - C;
err = 0;
for t=1:T,
    Z = exp(rho * log(Z) + shocks(t,:));
    W = (1 - delta) * K + Z.*K.^alpha;
    Cnew = x*W;
    err = err + abs(mean(beta * ((C ./ Cnew) .^ eta) .* (1 - delta + Z * alpha .* (K .^ (alpha - 1))) - 1));
    K = W - Cnew;
    C = Cnew;
end

function sc=Qscore(qq,C_t,W_t,K_t,Z_t,Z_next,Np)
if qq<0.001, sc=Qscore(0.001,C_t,W_t,K_t,Z_t,Z_next,Np); return, end
if qq>0.999, sc=Qscore(0.999,C_t,W_t,K_t,Z_t,Z_next,Np); return, end
global beta eta alpha delta N Cprev q
C_t(1:Np) = qq*W_t(1:Np);
eer_t_temp = beta * ((Cprev ./ C_t) .^ eta) .* (1 - delta + Z_t * alpha .* (K_t .^ (alpha - 1))) - 1;
eer_ts_temp = sum(eer_t_temp);

if eer_ts_temp<0,
    b=Np+ceil(rand*(N-Np));
    desired = eer_t_temp(b)-eer_ts_temp;
    C_t(b) = Cprev(b)*(((beta * (1 - delta + alpha*Z_t(b)*K_t(b)^(alpha-1)))/(1+desired))^(1.0/eta));
else
    sor=sortrows([eer_t_temp(Np+1:N);Np+1:N]',-1);
    i=1;
    while i<=10 && eer_ts_temp>1e-6,
        b=sor(i,2);
        desired = max(-1+1e-10,eer_t_temp(b)-eer_ts_temp);
        C_t(b) = Cprev(b)*(((beta * (1 - delta + alpha*Z_t(b)*K_t(b)^(alpha-1)))/(1+desired))^(1.0/eta));
        if C_t(b)>min(0.9,q+0.2)*W_t(b)
            C_t(b)=min(0.9,q+0.2)*W_t(b);
        end
        err_new_temp= beta * ((Cprev(b)/C_t(b))^eta) * (1-delta + Z_t(b)*alpha*(K_t(b)^(alpha - 1))) - 1;
        eer_ts_temp=eer_ts_temp-eer_t_temp(b)+err_new_temp;
        i=i+1;
    end
end

K_next = W_t - C_t;
W_next = (1-delta)*K_next+Z_next.*K_next.^alpha;
sc = sum(beta * ((C_t ./ W_next) .^ eta) .* (1 - delta + Z_next * alpha .* (K_next .^ (alpha - 1))) - 1);
