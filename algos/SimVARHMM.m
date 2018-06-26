function [x,reg] = SimVARHMM(mu,phi,sigma,Q,eta0,n)
%
% Generates a multivariate regime-switching autoregressive (p=1)
% Gaussian random walk
% regimes starting  from a given state eta0.
% 
% Input
%        mu    : d x reg matrix of means
%        phi   : d x reg autoregressive coefficient matrix
%        sigma : (d x d x reg)  sigma(:,:,k) is the volatility matrix 
%                                       of regime k
%        Q : transition matrix for the regimes;
%        eta0: initial state;
%        n   : length of the series.
%
% Output
%        x: series;
%        reg: simulated regimes;
%       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

d = size(mu,1);
x = zeros(n,d);

reg = SimMarkovChain(Q,n,eta0);

x(1,:) = mvnrnd(mu(:,reg(1)), sigma(:,:,reg(1)));

for i = 2:n
      x(i,:) = mvnrnd(mu(:,reg(i)), sigma(:,:,reg(i))) + (x(i-1,:) - mu(:,reg(i))') * phi(:,reg(i));
end

%%
function x = SimMarkovChain(Q,n,eta0)
%
% Generates a Markov chain X(1), ..., X(n) with transition matrix Q,
% starting from a state eta0 or the uniform distribution on {1,..., r}, 
% where r is the number of states.
%
%  Bruno Remillard, November 23, 2012.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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


