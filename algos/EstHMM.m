function [mu,sigma,Q,eta,nu,Z] = EstHMM(y,reg,max_iter,eps)
% Computes the estimation of parameters of a Gaussian univariate
% hidden Markov model.
%
%  INPUT
%   y: (n x 1) vector of returns;
%   reg: number of regimes;
%   prec: precision (stopping criteria); 
%               suggestion 0.0001 for daily returns;
%   max_iter: maxmimum number of iterations of the EM algo;
%
%
% OUTPUT
%   mu : 1 x reg vector of estimated means; 
%   sigma  : 1 x reg vector of estimated volatilities;
%   Q  : (reg x reg) estimated transition matrix;
%   eta: (n x reg) conditional probabilities of being in regime k at 
%        time t given observed returns r_1, ..., r_t;
%   nu : estimated stationary distribution.
%   Z : log likelihood of the joint distribution of the returns.
%
%
%    By Bruno Remillard, Nov 28, 2010
%    Improved by Christopher Leung, Sep 10, 2013
%%

ninit=100;  %minimum number of iterations
n = length(y);
r = reg;

%% Starting values for parameters
n0   = floor(n/r);
x    = zeros(n0,r);
ind0 = (1:n0)';

for j=1:r
    ind = (j-1)*n0 + ind0;
    x(:,j) = y(ind);
end

mu0    = mean(x);
sigma0 = std(x);
Q0     = ones(r,r)/r; 


%% warm-up
for k=1:ninit
  [nu, munew, sigmanew, Qnew, eta, gammabar, lambda, Lambda, Z] = EMStep(y,mu0,sigma0,Q0);
  mu0    = munew;
  Q0     = Qnew;
  sigma0 = sigmanew;
end

%% iterations
for k=1:max_iter
  [nu, munew, sigmanew, Qnew, eta, gammabar, lambda, Lambda, Z] = EMStep(y,mu0,sigma0,Q0);
  sum1 = sum(abs(mu0));
  sum2 = sum(abs(munew-mu0));
 
  if ( sum2 < sum1 * r * eps )
      break;
  end
  mu0    = munew;
  Q0     = Qnew;
  sigma0 = sigmanew;
  
end

%% output
mu    = munew;
sigma = sigmanew;
Q     = Qnew;

end

%%
function [nu, munew, sigmanew, Qnew, eta, gammabar, lambda, Lambda, Z] = EMStep(y,mu,sigma,Q)

n = length(y);
r = length(mu);
%gamma    = zeros(n,r);
gammabar = zeros(n,r);
lambda   = zeros(n,r);
f        = zeros(n,r);
f_       = zeros(n,r);
Lambda   = zeros(r,r,n);
M = zeros(r,r);

%% f
for j=1:r
  z      = (y-mu(j))/sigma(j);
  f_(:,j) = normpdf(z);
  f(:,j) = f_(:,j) / sigma(j);
end

%%
gammabar(n,:)=1/r;

for k = 1: (n-1)
    i = n-k;
    j = i+1;
    v =  ( gammabar(j,:) .* f(j,:) ) * Q' ;
    gammabar(i,:) = v/sum(v); %normalized so that it does not explode
end

%%
eta0 = ones(1,r)/r;
eta = zeros(n,r);
Z_cond = zeros(n,1);

v = ( eta0 * Q) .* f(1,:);
Z_cond(1) = sum(v);
eta(1,:) = v/Z_cond(1);

for i=2:n
    v        = ( eta(i-1,:) * Q) .* f(i,:);
    Z_cond(i) = sum(v);
    eta(i,:) = v/Z_cond(i); 
end

Z = sum(log(Z_cond));

%% 10.14 (in Bruno's book)
v      = eta .* gammabar ;
sv0    = sum(v,2);
% sv     = repmat( sv0 , 1, r);
% lambda = v ./ sv ;

for j=1:r
    lambda(:,j) = v(:,j) ./ sv0;
end

%% 10.15 (in Bruno's book)

 for j=1:r
     Lambda(j,:,n) = lambda(n,j) * Q(j,:) ;
 end
 
gf =  gammabar .* f ;

for i=1:(n-1)
    M = Q .* ( eta(i,:)' * gf(i+1,:)  );
    c = sum(sum(M)) ;
    Lambda(:,:,i) = M/c;
end

%% M-Step
nu = mean(lambda);
munew = zeros(1,r);
sigmanew = zeros(1,r);
Qnew = zeros(r,r);

for j=1:r 
         w  = lambda(:,j) / nu(j) /n;
   munew(j) = sum( y .* w );
sigmanew(j) = sqrt( sum( y.^2 .* w) - munew(j)^2 );
  Qnew(j,:) = mean(Lambda(j,:,:),3) / nu(j) ;
end

end