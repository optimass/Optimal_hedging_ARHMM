function out = GofVHMM(R,reg,max_iter, prec,n_sample)
% Goodness-of-fit test for multivariate Gaussian HMM models 
% using Rosenblatt's transform and the parametric bootstrap technique
%
%  Input: R (n x p) matrix of returns
%         n: lenght of series, p: dimension
%         reg: number of regimes
%         prec: precision (stopping criteria)
%               suggestion 0.0001 for daily returns;
%         max_iter: maxmimum number of iterations of the EM algo
%         n_sample: number of simulated samples for parametric bootstrap
%
%
%  output: mu : p x reg matrix of estimated means 
%          A  : (p x p x reg)  A(:,:,k) is the covariance matrix 
%                                       of regime k
%          Q  : (reg x reg) estimated transition matrix
%          nu : estimated stationary distribution
%          eta: (n x reg) eta(t,k) is the conditional 
%                         prob of being in regime k at time t
%                         given observed returns r_1, ..., r_t
%          cvm_est: Cramer-von Mises statistic for the data
%          cvm_sim: Cramer-von Mises statistics for the simulated samples
%          pvalue: p-value (in percent) for the Cramer-von Mises statistic.
%
%     By Bruno Remillard and Massimo Caccia, Apr 3rd, 2015
%%

[mu,A,Q,nu,eta,cvm_est] = EstVHMM(R,reg,max_iter, prec);

fprintf('\n End of estimation\n'); 

[n,~] = size(R) ; % #of obs and dimensions


%% Parametric bootstrap

cvm_sim = zeros(n_sample,1);
eta0 = randsample(reg,n_sample,true);

tic
parfor i=1:n_sample
    
    R1 = SimVHMM(mu,A,Q,eta0(i),n);
    
    [~,~,~,~,~,cvm_sim(i)] = EstVHMM(R1,reg,max_iter, prec);
    
end

toc

%% output

pvalue = 100*mean( cvm_sim > cvm_est);

out.pvalue  = pvalue;
out.mu      = mu;
out.A       = A;
out.Q       = Q;
out.nu      = nu;
out.eta     = eta;
out.cvm_est = cvm_est;
out.cvm_sim = cvm_sim;


