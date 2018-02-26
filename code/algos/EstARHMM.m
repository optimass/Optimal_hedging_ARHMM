function [mu,phi,sigma,Q,nu,eta,Z] = EstARHMM(y,r,max_iter,prec) %codegen
% Computes the estimation of parameters of a Gaussian univariate
% autoregressive (p=1) regime-switching model.
%
%  INPUT
%   y: (n x 1) vector of returns;
%   r: number of regimes;
%   max_iter: maxmimum number of iterations of the EM algo;
%   prec: precision (stopping criteria); suggestion 0.0001 for daily returns;
%
%
% OUTPUT
%   mu    : reg x 1 vector of estimated means; 
%   phi   : reg x 1 vector of estimated autoregressive coefficient;
%   sigma : reg x 1 vector of estimated volatilities;
%   Q     : (reg x reg) estimated transition matrix;
%   eta: (n x reg) conditional probabilities of being in regime k at 
%        time t given observed returns r_1, ..., r_t;
%   nu    : estimated stationary distribution.
%   Z     : log likelihood of the joint distribution of the returns.
%
%
%    By Bruno Remillard and Massimo Caccia, Apr 3rd, 2015
%%

ninit=100;  %minimum number of iterations
n = size(y,1); %number of obs


%% Starting values for parameters

n0   = floor(n/r);
x    = zeros(n0,r);
ind0 = (1:n0)';

for j=1:r
    ind = (j-1)*n0 + ind0;
    x(:,j) = y(ind);
end

mu0 = mean(x)';
sigma0 = std(x)';
phi0 = zeros(r,1);
Q0 = ones(r,r)/r; 

%% Part 1 :: ar0

ar1 = 0;

% warm-up

for k=1:ninit
  [nu, munew, phinew, sigmanew, Qnew, eta, ~, ~, ~, Z] = EMStep(y,mu0,phi0,sigma0,Q0,ar1);
  mu0    = munew;
  Q0     = Qnew;
  sigma0 = sigmanew;
  phi0   = phinew;
end

%  iterations

for k=1:max_iter
  [nu, munew, phinew, sigmanew, Qnew, eta, ~, ~, ~, Z] = EMStep(y,mu0,phi0,sigma0,Q0,ar1);
  sum1 = sum(abs(mu0));
  sum2 = sum(abs(munew-mu0));
 
  if ( sum2 < sum1 * r * prec )
      break;
  end
  
  mu0    = munew;
  Q0     = Qnew;
  sigma0 = sigmanew;
  phi0   = phinew;
end

%% Part 2 :: ar1

ar1 = 1;

%  iterations

for k=1:max_iter
  [nu, munew, phinew, sigmanew, Qnew, eta, ~, ~, ~, Z] = EMStep(y,mu0,phi0,sigma0,Q0,ar1);
  sum1 = sum(abs(phi0));
  sum2 = sum(abs(phinew-phi0));
 
  if ( sum2 < sum1 * r * prec )
      break;
  end
  
  mu0    = munew;
  Q0     = Qnew;
  sigma0 = sigmanew;
  phi0   = phinew;
end


%% output

mu  = munew;
sigma = sigmanew;
phi = phinew;  
Q   = Qnew;

end


function [nu, munew, phinew, sigmanew, Qnew, eta, gammabar, lambda, Lambda, Z] = EMStep(y,mu,phi,sigma,Q,ar1)

n = length(y)-1; %we loose one obs because of the AR(1)
r = size(mu,1);

gammabar = zeros(n,r);
lambda   = zeros(n,r); % 10.14
f        = zeros(n,r);
Lambda   = zeros(r,r,n); % 10.15

%%

for j=1:r
    z = ( y(2:end) - (mu(j)+ phi(j)*(y(1:end-1)-mu(j)) ) ) / sigma(j);
    f(:,j) = normpdf(z)/ sigma(j);

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

%%
v      = eta .* gammabar ;
sv0    = sum(v,2);
% sv     = repmat( sv0 , 1, r);
% lambda = v ./ sv ;

for j=1:r
    lambda(:,j) = v(:,j) ./ sv0;
end

%%

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

phinew = zeros(r,1);  munew = zeros(r,1);
sigmanew = zeros(r,1);    Qnew = zeros(r); 


for j=1:r
    
    w  = lambda(:,j) / nu(j) /n;
    ybar = y(2:end)' * w;
    y_ = y(1:end-1)' * w;
    
    
    if ar1
        temp1 = ( (y(1:end-1)-y_).*sqrt(w) )' * ( (y(1:end-1)-y_).*sqrt(w) );
        temp2 = ( (y(2:end)-ybar).*sqrt(w) )' * ( (y(1:end-1)-y_).*sqrt(w) );
        phinew(j) = temp2 / temp1; % (15)
    end

    munew(j) = (ybar-phinew(j)*y_) / (1-phinew(j));  % (14)
    
    e = y(2:end) - ybar - phinew(j) * ( y(1:end-1) - y_ );
    
    sigmanew(j) = sqrt((e.*sqrt(w))' * (e.*sqrt(w))); % (16)
    
    Qnew(j,:) = mean(Lambda(j,:,:),3) / nu(j); % (13)

end


end

