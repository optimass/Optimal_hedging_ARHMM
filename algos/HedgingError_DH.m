function [HE_BS,C0] = HedgingError_DH(S,T,K,r,mu,sigma,put)
%
%   Input:
%       mu : annualized
%       sigma: anunalized
%
%
%
%
S0 = S(1);
n = length(S)-1;

% for simulation: no dividend
q=0;

Tp = T/n;
rp = r*Tp;
Kp = K*exp(-r*T);
F = S0*exp(-q*T);


[N1, N2] = computed1d2(S0,K,mu,q,sigma,T,0,put);

if put
    C0 = (Kp*N2 - F*N1);
else
    C0 = (F*N1 - Kp*N2);
end


%% need to start the replication

% hedging portfolio
V = C0;

Stilde0 = S0;
t = 0;

for j = 1:n
    Stilde1 = exp(-rp*j)*S(j+1);
    Delta =   Stilde1-Stilde0;
    phiS =    computed1d2(S(j),K,mu,q,sigma,T,t,put);    
    V = V+phiS*Delta;
    Stilde0 = Stilde1;
    t = t+Tp;      
end

if put
  payoff = max(0,Kp-Stilde1);
else
  payoff = max(0,Stilde1-Kp);
end


HE_BS =  payoff-V;


function [N1, N2] = computed1d2(S0,K,mu,q,sigma,T,t,put)

d1 = log(S0/K) + (mu-q+0.5*sigma^2)*(T-t)  ;
d1 = d1 / ( sigma*sqrt(T-t) );
d2 = d1 - sigma*sqrt(T-t);


if put
    N1 = normcdf(-d1);
    N2 = normcdf(-d2);
else
    N1 = normcdf(d1);
    N2 = normcdf(d2);
end
















