function out =  GofARHMM(R,reg,max_iter,prec,N)
% Goodness-of-fit test for autoregressive Gaussian HMM models 
% using Rosenblatt's transform and the parametric bootstrap technique
%
%  Input: R (n x 1) matrix of returns
%         n: lenght of series, p: dimension
%         reg: number of regimes
%         prec: precision (stopping criteria)
%               suggestion 0.0001 for daily returns;
%         max_iter: maxmimum number of iterations of the EM algo
%         N: number of simulated samples for parametric bootstrap
%
%
%  output: mu : reg x 1 vector of estimated means 
%          A  : reg x 1 vector of estimated variance
%          phi: reg x 1 vector of estimated autoregressive parameters
%          Q  : (reg x reg) estimated transition matrix
%          nu : estimated stationary distribution
%          eta: (n x reg) eta(t,k) is the conditional 
%                         prob of being in regime k at time t
%                         given observed returns r_1, ..., r_t
%          cvm_est: Cramer-von Mises statistic for the data
%          cvm_sim: Cramer-von Mises statistics for the simulated samples
%          pvalue: p-value (in percent) for the Cramer-von Mises statistic.
%
%
%     By Bruno Remillard and Massimo Caccia Apr 15, 2015
%%

n = length(R);
Z = zeros(n-1,reg);

[mu,phi,sigma,Q,nu,eta] = EstARHMM(R,reg,max_iter,eps);

out.mu    = mu;
out.sigma = sigma;
out.phi   = phi;
out.Q     = Q;
out.eta   = eta;
out.nu    = nu;

for j=1:reg
    Z(:,j) = ( R(2:end) - (mu(j)+ phi(j)*(R(1:end-1)-repmat(mu(j),n-1,1))) ) / sigma(j) ; 
end

U = normcdf(Z);

eta00 = ones(1,reg)/reg;

w00 = [eta00;eta] * Q; % ??
w   = w00(1:end-1,:);  % le w du livre
W = sum( w .* U,2);    % le PSI du livre

cvm = SnB(W);

cvm_sim = zeros(N,1);

parfor i=1:N

    R1 = SimARHMM(mu,phi,sigma,Q,reg,n);
    
    [mu1,phi1,sigma1,Q1,~,eta1] = EstARHMM(R1,reg,max_iter,prec);
 
    Z1 = (repmat(R1(2:end),1,reg) - (repmat(mu1',n-1,1) + repmat(phi1',n-1,1).*( repmat(R1(1:end-1),1,reg) - (repmat(mu1',n-1,1) ))))...
                    ./ repmat(sigma1',n-1,1); 

    U = normcdf(Z1);


    w00 = [eta00;eta1] * Q1;
    w   = w00(1:end-1,:);
    W1 = sum( w .* U,2);

cvm_sim(i) = SnB(W1);
       
end

pvalue = 100*mean(cvm_sim>cvm);

out.cvm = cvm;
out.pvalue = pvalue;




function  stat = SnB(U)

n = length(U);
l = (2*(1:n)-1)/n;

U0 = sort(U);

stat = n/3 +U0'*U0 - l*U0;  % comprend pas



function [x,reg] = SimHMMVAR1d(mu,phi,sigma,Q,eta0,n)
%
% Generates a regime-switching random walk with autoregressive (p=1) Gaussian
% regimes starting  from a given state eta0.
% 
% Input
%        mu: mean of each regime;
%        A : variance of each regime;
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


