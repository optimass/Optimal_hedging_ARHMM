function out =  GofHMM1d(R,reg,max_iter,prec,N)


n = length(R);
Z = zeros(n,reg);

[mu,sigma,Q,eta,nu] = EstHMM1d_mex(R,reg,max_iter,prec);

out.mu    = mu;
out.sigma = sigma;
out.Q     = Q;
out.eta   = eta;
out.nu    = nu;

for j=1:reg
Z(:,j) = (R-mu(j) ) / sigma(j) ; 
end

U = normcdf(Z);

eta00 = ones(1,reg)/reg;

w00 = [eta00;eta] * Q;

w   = w00(1:end-1,:);

W = sum( w .* U,2);

cvm = SnB(W);

cvm_sim = zeros(N,1);

parfor i=1:N
    
    R1 = SimHMMGaussian1d(mu,sigma,Q,reg,n);
    
    
   [mu1,sigma1,Q1,eta1,nu1] = EstHMM1d_mex(R1,reg,max_iter,prec);


Z1 = (repmat(R1,1,reg) - repmat(mu1,n,1)) ./ repmat(sigma1,n,1); 

U = normcdf(Z1);


w00 = [eta00;eta1] * Q1;

w   = w00(1:end-1,:);

W1 = sum( w .* U,2);

cvm_sim(i,1) = SnB(W1);
    
end


cvm
pvalue = 100*mean(cvm_sim>cvm)

out.cvm = cvm;
out.pvalue = pvalue;

end
%%
function  stat = SnB(U)

n = length(U);
l = (2*(1:n)-1)/n;

U0 = sort(U);

stat = n/3 +U0'*U0 - l*U0;  % CvM stat
end

%%
function [x,reg] = SimHMMGaussian1d(mu,sigma,Q,eta0, n)
%
% Generates a regime-switching random walk with Gaussian regimes starting 
% from a given state eta0.
% 
% Input
%        mu: mean of each regime;
%        sigma : volatilies of each regime;
%        Q : transition matric for the regimes;
%        eta0: initial state;
%        n   : length of the series.
%
%      Output
%         x: series;
%        reg: simulated regimes;
%       
%
%  Bruno Remillard, November 23, 2012.

r = length(mu);
x = zeros(n,1);



reg = SimMarkovChain(Q,n,eta0);

  x0 = zeros(n,r);
    
  for k=1:r
     x0(:,k) = mu(k)+ sigma(k)*randn(n,1);
  end
  
  for i=1:n
      x(i) = x0(i,reg(i) );
  end

end
%%
function x = SimMarkovChain(Q,n,eta0)
%
% Generates a Markov chain X(1), ..., X(n) with transition matrix Q,
% starting from a state eta0 or the uniform distribution on {1,..., r}, 
% where r is the number of states.
%
%  Bruno Remillard, November 23, 2012.



[r,p] = size(Q);
x = zeros(n,1);
x0 = zeros(n,r);

if(nargin <3)
    ind = randsample(r,1);
else
    ind = eta0;
end

for k=1:r
   x0(:,k) = randsample(r,n,true, Q(k,:) );
end

for i=1:n
   x(i) =  x0(i,ind);
   ind = x(i);
end
end