function [mu,A,Q,eta,nu,cvm] = EstVHMM(y,r,max_iter,prec)
% Computes the estimation of parameters of a Gaussian multivariate
% regime-switching model.
%
%  INPUT
%   y: (n x d) vector of returns;
%   r: number of regimes;
%   max_iter: maxmimum number of iterations of the EM algo;
%   prec: precision (stopping criteria); suggestion 0.0001 for daily returns;
%
%
% OUTPUT
%   mu    : d x reg matrix of estimated means
%   A     : (d x d x reg)  sigma(:,:,k) is the covariance matrix 
%                                       of regime k
%   Q     : (reg x reg) estimated transition matrix;
%   eta: (n x reg) conditional probabilities of being in regime k at 
%        time t given observed returns r_1, ..., r_t;
%   nu : estimated stationary distribution.
%   cvm: Cramer-von-Mises statistic for goodness-of-fit.
%
%    By Bruno Remillard and Massimo Caccia, Apr 3rd, 2015
%%

ninit=100;  %minimum number of iterations
n = size(y,1); %number of obs
d = size(y,2); %number of variables


%% Starting values for parameters

n0   = floor(n/r);      ind0 = (1:n0)';
x    = zeros(n0,r,d);   X = zeros(n0,d,r);
mu0  = zeros(d,r);      sigma0 = zeros(d,d,r);

for i=1:d
    for j=1:r
        ind = (j-1)*n0 + ind0;
        x(:,j,i) = y(ind,i); %to compute mu
        X(:,i,j) = y(ind,i); %to compute sigma
    end
end

for i=1:d
    mu0(i,:) = mean(x(:,:,i));
end

for i=1:r
    sigma0(:,:,i) = cov(X(:,:,i));
end

Q0 = ones(r,r)/r;

%% warm-up

for k=1:ninit
  [nu, munew, sigmanew, Qnew, eta, ~, ~, ~] = EMStep(y,mu0,sigma0,Q0);
  mu0    = munew;
  Q0     = Qnew;
  sigma0 = sigmanew;
  
end

%% iterations

for k=1:max_iter
  [nu, munew, sigmanew, Qnew, eta, ~, ~, ~] = EMStep(y,mu0,sigma0,Q0);
  sum1 = sum(sum(abs(mu0)));
  sum2 = sum(sum(abs(munew-mu0)));

  if ( sum2 < sum1 * r * prec )
      break;
  end
  mu0    = munew;
  Q0     = Qnew;
  sigma0 = sigmanew;
  
end
%% CVM test

Psi = RosenblattHMM(y,mu0,sigma0,r,n,d,eta,Q0);

cvm = SnB(Psi);

%% output

mu = munew;
A  = sigmanew;
Q  = Qnew;

end

%%
function [nu, munew, sigmanew, Qnew, eta, gammabar, lambda, Lambda] = EMStep(y,mu,sigma,Q)

n = length(y);
r = size(mu,2);
d = size(mu,1);
gammabar = zeros(n,r);
lambda   = zeros(n,r); % 10.14
f        = zeros(n,r);
Lambda   = zeros(r,r,n); % 10.15


%%

for j=1:r
  f(:,j) = vnpdf(y,mu(:,j)',sigma(:,:,j));
  %f(:,j) = normpdf(:,z)/sigma(:,j);
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
eta0 = ones(1,r)/r; eta = zeros(n,r);

v = ( eta0 * Q) .* f(1,:);
eta(1,:) = v/sum(v);

for i=2:n
    v        = ( eta(i-1,:) * Q) .* f(i,:);
    eta(i,:) = v/sum(v);
end

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

munew = zeros(d,r); sigmanew = zeros(d,d,r) ; Qnew = zeros(r,r);

for j=1:r 
w  = lambda(:,j) / nu(j) /n;
munew(:,j) = sum( y .* (w*ones(1,d)) );

z = y - ones(n,1)*munew(:,j)';
sigmanew(:,:,j) = (z.*(sqrt(w)*ones(1,d)))' * (z.*(sqrt(w)*ones(1,d)));

Qnew(j,:) = mean(Lambda(j,:,:),3) / nu(j) ;

end

%%%% singularity problem %%%%
eps=100;
for j=1:r
    eps = min( min(diag(sigmanew(:,:,j)), eps) );
end
sigmanew = sigmanew + eps*0.0001;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

function out = vnpdf(X,mu,A)

n = length(X);
out = zeros(n,1);
z = X - repmat(mu,n,1);
invA = inv(A);
detA = det(A);

for i = 1:n
    out(i,1) = exp( - z(i,:) * invA * z(i,:)' / 2 ) / sqrt(detA) ;
end

end




%%
function Psi = RosenblattHMM(R,mu,A,reg,n,p,eta,Q)
% function to compute Rosenblatt's Transform

Z = zeros(n,p,reg);
for i = 1:reg
    Z(:,1,i) = (R(:,1) - repmat(mu(1,i),n,1)) / sqrt(A(1,1,i)) ; % first dimension
    
    for j = 2:p % Rosenblatt Gaussian

       R11 = A(1:j-1,1:j-1,i);
       R21 = A(j,1:j-1,i);
       R11inv = inv(R11);
       B = R21*R11inv;
       
       x = R(:,1:j-1);
       y = R(:,j);
       
       Omega = A(j,j,i) - B*R11*B'; 
       m = repmat( mu(j,i)',n,1 ) + ( x - repmat( mu(1:j-1,i)',n,1) ) * B' ; 
   
       Z(:,j,i) = (y-m) / sqrt(Omega);
       
    end 
end 

F = normcdf(Z);

%%

eta00 = ones(1,reg)/reg;  
w00 = [eta00;eta] * Q;    
W   = w00(1:end-1,:);    

%% 

Psi = zeros(n,p);   % prealloc
temp = zeros(n,reg);% prealloc
temp2 = temp ;      % prealloc

for r = 1:reg
    temp(:,r) = F(:,1,r) ;
end
Psi(:,1) =  sum( W .* temp,2);

for q = 2:p      
    for r = 1:reg
        %temp(:,r)  = W(:,r) .* mvnpdf( R(:,1:q-1) , mu(1:q-1,r)' , A(1:q-1,1:q-1,r) ) ;
        % we need vnpdf bc of the matlab coder:
        temp(:,r)  = W(:,r) .* vnpdf( R(:,1:q-1) , mu(1:q-1,r)' , A(1:q-1,1:q-1,r) ) ;
        temp2(:,r) = F(:,q,r);
    end
    Psi(:,q) = sum( temp .* temp2 , 2)  ./ sum(temp,2) ;
end

end
%%
function  Sn = SnB(E)
%  Cramer-von Mises statistic SnB for GOF based on the Rosenblatt transform
%  Ref: Genest, Remillard & Beaudoin 2009
%
%   Input
%        E: n x d matrix of pseudo-observations.
%
%   Output
%        Sn: Cramer-von Mises statistic SnB.

[n,d] = size(E);

Dn = zeros(n,1);
S1 = n/3^d;

G0 = 1-E.*E;
E0 = 1-E;

S2 = sum(prod(G0,2))/2^(d-1);


for i=1:n
    G0 = repmat(E0(i,:),n,1);
    
    Dn(i)= mean( prod( min(G0,E0),2));
end

Sn = sum(Dn)-S2+S1;

end


