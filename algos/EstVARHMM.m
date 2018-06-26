function [mu,phi,A,Q,nu,eta,cvm] = EstVARHMM(y,r,max_iter,eps) %#codegen
% Computes the estimation of parameters of a Gaussian multivariate
% autoregressive ( p=1 ) regime-switching model.
%
%  INPUT
%   y: (n x d) vector of returns;
%   r: number of regimes;
%   prec: precision (stopping criteria); suggestion 0.0001 for daily returns;
%   max_iter: maxmimum number of iterations of the EM algo;
%
%
%  OUTPUT
%          mu : p x reg matrix of estimated means 
%          phi: (p x p x reg)  PHI(:,:,k) is the autoregressive coefficent
%                                       of regime k
%          A  : (p x p x reg)  A(:,:,k) is the covariance matrix 
%                                       of regime k
%          Q  : (reg x reg) estimated transition matrix
%          nu : estimated stationary distribution
%          eta: (n x reg) conditional probabilities of being in regime k at 
%                           time t given observed returns r_1, ..., r_t;
%          nu : estimated stationary distribution.
%          cvm: Cramer-von-Mises statistic for goodness-of-fit.
%
%    By Bruno Remillard and Massimo Caccia, Apr 3, 2015
%%

ninit=100;     %minimum number of iterations
[n,d] = size(y); %number of obs and dimensions


%% Starting values for parameters

n0 = floor(n/r);      ind0 = (1:n0)';

x   = zeros(n0,r,d);  X = zeros(n0,d,r);      %prealloc
mu0 = zeros(d,r);     sigma0 = zeros(d,d,r);      phi0 = zeros(d,d,r);

for i=1:d
    for j=1:r
        ind = (j-1)*n0 + ind0;
        x(:,j,i) = y(ind,i); %to compute mu
        X(:,i,j) = y(ind,i); %to compute A
    end
end

for i = 1:d
    mu0(i,:) = mean(x(:,:,i));
end

for i = 1:r
    sigma0(:,:,i)   = cov( X(:,:,i) );
    phi0(:,:,i) = X(1:end-1,:,i) \ X(2:end,:,i) ; 
end

Q0 = ones(r,r)/r;

%% warm-up

for k = 1:ninit
  [nu, munew, sigmanew, phinew, Qnew, eta, ~, ~, ~] = EMStep(y,mu0,sigma0,phi0,Q0);
  mu0    = munew;
  Q0     = Qnew;
  sigma0 = sigmanew;
  phi0   = phinew;
end

%% iterations

for k = 1:max_iter
  [nu, munew, sigmanew, phinew, Qnew, eta, ~, ~, ~] = EMStep(y,mu0,sigma0,phi0,Q0);
  sum1 = sum(sum(abs(mu0)));
  sum2 = sum(sum(abs(munew-mu0))); 

  if ( sum2 < sum1 * r * eps * d ) % stopping criteria
      break;
  end
  
  mu0    = munew;
  Q0     = Qnew;
  sigma0 = sigmanew;
  phi0   = phinew;
end

%% CVM test

Psi = RosenblattHMMVAR(y,mu0,sigma0,phi0,r,n,d,eta,Q0);

cvm = SnB(Psi);

%% output

mu    = munew;
A     = sigmanew;
phi   = phinew;  
Q     = Qnew;

end




function [nu, munew, sigmanew, phinew, Qnew, eta, gammabar, lambda, Lambda] = EMStep(y,mu,sigma,phi,Q)

n        = length(y)-1; %we loose one obs because of the AR(1)
[d,r]    = size(mu);
% gamma  = zeros(n,r);
gammabar = zeros(n,r);
lambda   = zeros(n,r); % 10.14
f        = zeros(n,r);
Lambda   = zeros(r,r,n); % 10.15

%% pdf

for j=1:r
    f(:,j) = varpdf(y,mu(:,j)',sigma(:,:,j),phi(:,:,j));
   %f(:,j) = normpdf(:,z)/sigma(:,j);
end

%% gammabar

gammabar(n,:) = 1/r;

for k = 1: (n-1)
    i = n-k;
    j = i+1;
    v = ( gammabar(j,:) .* f(j,:) ) * Q' ;
    gammabar(i,:) = v/sum(v); %normalized so that it does not explode
end


%% eta

eta0 = ones(1,r)/r;   eta = zeros(n,r) ;

v = ( eta0 * Q) .* f(1,:);
eta(1,:) = v/sum(v);

for i=2:n
    v        = ( eta(i-1,:) * Q) .* f(i,:);
    eta(i,:) = v/sum(v);
end

%% lambda

v      = eta .* gammabar ;
sv0    = sum(v,2);
% sv     = repmat( sv0 , 1, r);
% lambda = v ./ sv ;

for j=1:r
    lambda(:,j) = v(:,j) ./ sv0;
end

%% Lambda

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

munew = zeros(d,r);        phinew = zeros(d,d,r);    %prealloc
sigmanew = zeros(d,d,r);   Qnew = zeros(r);

for j=1:r
    
    w    = lambda(:,j) / nu(j) /n;
    ybar = y(2:end,:)'   * w;
    y_   = y(1:end-1,:)' * w;

    temp1 = ( (y(1:end-1,:)- repmat(y_',n,1))   .* (sqrt(w)*ones(1,d)) )' * ...
                            ( (y(1:end-1,:)-repmat(y_',n,1)) .* (sqrt(w)*ones(1,d)) );
    temp2 = ( (y(2:end,:)  - repmat(ybar',n,1)) .* (sqrt(w)*ones(1,d)) )' * ...
                            ( (y(1:end-1,:)-repmat(y_',n,1)) .* (sqrt(w)*ones(1,d)) );
    
                        
    phinew(:,:,j) = temp2*pinv(temp1); % (15)
    
    munew(:,j) = pinv( eye(d) - phinew(:,:,j) ) * (ybar - phinew(:,:,j)*y_); % (14)
    
    e = y(2:end,:) - repmat(ybar',n,1) - ( y(1:end-1,:) - repmat(y_',n,1) ) * phinew(:,:,j)';
    sigmanew(:,:,j) = (e.*(sqrt(w)*ones(1,d)))' * (e.*(sqrt(w)*ones(1,d))); % (16)
    
    Qnew(j,:) = mean(Lambda(j,:,:),3) / nu(j) ; % (13)

end

%%%% singularity problem %%%%
eps=100;
for j=1:r
    eps = min( min(diag(sigmanew(:,:,j)), eps) );
end
sigmanew = sigmanew + eps*0.0001;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end


function Psi = RosenblattHMMVAR(R,mu,A,phi,reg,n,p,eta,Q)
% Computes the Rosenblat transform for HMMVAR
%
%% create a new MU

n = n-1; % we lose one obs
MU = zeros(n,p,reg); 

for r = 1:reg 
    for i = 1:n
        MU(i,:,r) = mu(:,r)' + (R(i,:)-mu(:,r)') * phi(:,:,r)' ;
    end
end    

R(1,:) = [] ;

%% create F

Z = zeros(n,p,reg);

for i = 1:reg
    Z(:,1,i) = (R(:,1) - MU(:,1,i) ) / sqrt(A(1,1,i)) ; % first dimension
    
    for j = 2:p % Rosenblatt Gaussian

       R11 = A(1:j-1,1:j-1,i);
       R21 = A(j,1:j-1,i);
       R11inv = inv(R11);
       B = R21*R11inv;
       
       x = R(:,1:j-1);
       y = R(:,j);
       
       Omega = A(j,j,i) - B*R11*B';
       m = MU(:,j,i) + ( x - MU(:,1:j-1,i) ) * B' ; 
   
       Z(:,j,i) = (y-m) / sqrt(Omega);
       
    end 
end 

F = normcdf(Z);

%% create W

eta00 = ones(1,reg)/reg; % 
w00 = [eta00;eta] * Q;   % 
W   = w00(1:end-1,:);    % 

%% Create Psi 

Psi = zeros(n,p);   temp = zeros(n,reg);   temp2 = temp ;    % prealloc

for r = 1:reg
    temp(:,r) = F(:,1,r) ;
end
Psi(:,1) =  sum( W .* temp , 2);

for q = 2:p         
    for r = 1:reg
        temp(:,r)  = W(:,r) .* vnpdf( R(:,1:q-1) , MU(:,1:q-1,r) , A(1:q-1,1:q-1,r) ) ;
        temp2(:,r) = F(:,q,r);
    end
    Psi(:,q) = sum( temp .* temp2 , 2)  ./ sum(temp,2) ;
end

end


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



function y = varpdf(X,MU,A,PHI)
% Vector Autoregressive (p=1) probability density function (proportional)
%
n = size(X,1)-1; % one obs is lost

z = X(2:end,:) - repmat(MU,n,1) - ( X(1:end-1,:) - repmat(MU,n,1) ) * PHI';

Ainv = pinv(A);

y = zeros(n,1);

for i = 1:n
    y(i) = exp( -0.5 * z(i,:) * Ainv * z(i,:)');
end

y = y  / sqrt(det(A)) ;
    
end


function out = vnpdf(X,mu,A)
% Vector gaussian probability density function (proportional!)
% SPECIAL because mu is nxd

n    = length(X);
out  = zeros(n,1);
z    = X - mu;
invA = inv(A);
detA = det(A);

for i = 1:n
    out(i,1) = exp( - z(i,:) * invA * z(i,:)' / 2 ) / sqrt(detA) ;
end

end

