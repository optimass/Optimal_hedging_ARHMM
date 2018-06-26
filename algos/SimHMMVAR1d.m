function [x,reg] = SimHMMVAR1d(mu,phi,sigma,Q,eta0,n)
%
% Generates a regime-switching random walk with autoregressive (p=1) Gaussian
% regimes starting  from a given state eta0.
% 
% Input
%        mu: mean of each regime;
%        sigma : vol of each regime;
%        Q : transition matric for the regimes;
%        eta0: initial state;
%        n   : length of the series.
%
% Output
%        x: series;
%        reg: simulated regimes;
%       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r = size(mu,1);
x = zeros(n,1);

reg = SimMarkovChain(Q,n,eta0);

x(1) = mu(reg(1)) +  randn(1) * sigma(reg(1));

for i = 2:n
    x(i) = mu(reg(i)) +  randn(1) * sigma(reg(i)) + (x(i-1) - mu(reg(i))) * phi(reg(i));
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


