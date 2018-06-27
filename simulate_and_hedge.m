%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% In this script, we will simulate processes and hedge the options 
% according to the following mthods: 
%
% 1) Black-Scholes (BS)
% 
% 2) Optimal Hedging assuming returns are Gaussian (OH-BS)
%
% 3) Optimal Hedging assuming returns follow an regime switching random
%    walk (OH-HMM)
%
% 4) Optimal Hedding assuming returns follow and autoregressive regime
%    switching random walk (OH-ARHMM)
%
% models from the paper :
%       "Option Pricing and Hedging for Discrete 
%        Time Autoregressive Hidden Markov Model"
%
% https://papers.ssrn.com/sol3/papers.cfm?abstract_id=2995944
%
% This script reproduces Figure 8 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parameters

reg = 3;                  % number of regimes
max_iter = 1000;          % max number of iterations when calibrating
prec = 0.0001;            % precision threshold when calibrating
N = 10000;                % number of samples for Monte Carlo
n = 63;                   % duration of the option (in days)
T = n/252;                % duration (in years)
f = 0;                    % forward rate
r = 0.01;                 % interest rate
train_win = 10000;        % training window for calibration

S0 = 100;                 % starting stock price
K = 100;                  % strike price
put=0;                    % put or call 
method = 1 ;              % use "most probable regime" when hedging

m = 100;                  % number of points on the grid
q = 30;                   % q: number of points on the lag return grid
minS = 65;                % min Stock price on the grid
maxS = 135;               % max Stock price on the grid



%% Set the ARHMM parameters:

mu_annual = [0.3 0.1 -0.3];
phi = [0.2 0 -0.2];
sigma_annual = [0.2 0.4 0.8];
Q = [0.95  0.025 0.025 ; 
     0.025 0.95  0.025 ;
     0.025 0.025 0.95];
eta0 = 1;

% de-anualize:
mu = mu_annual / 252 ;
sigma = sigma_annual / sqrt(252);

%% Simulate a process to calibrate a model for each Methodologies:

R_train = SimARHMM(mu,phi,sigma,Q,eta0,train_win);


%% estimate params for hedging

[mu_hmmvar,phi_hmmvar,sigma_hmmvar,Q_hmmvar,nu_hmmvar,eta_hmmvar] = EstARHMM(R_train,reg,max_iter,eps);
[mu_hmm,sigma_hmm,Q_hmm,nu_hmm,eta_hmm]                           = EstHMM(R_train,reg,max_iter,eps);
[mu_bs,sigma_bs,~,~,~]                                            = EstHMM(R_train,1,max_iter,eps);

% anualize parameters:
mu_bs     = 252*mu_bs        ; sigma_bs     = sqrt(252)*sigma_bs     ;
mu_hmm    = 252*mu_hmm'      ; sigma_hmm    = sqrt(252)*sigma_hmm'   ;
mu_hmmvar = 252*mu_hmmvar    ; sigma_hmmvar = sqrt(252)*sigma_hmmvar ;



%% align the starting regime for simulation and for hedging

eta0_simu = 1;
[~, idx] = max(mu_hmmvar); % always the same as hmm
eta0 = zeros(reg,1);
eta0(idx) = 1;


%% Create hedging grids

% OH-BS
[S1,C1,a1,gamma1,rho1] = HedgingGaussian(mu_bs,sigma_bs,T,K,r,n,put,minS,maxS,S0,m);
display('OH-BS done');

% OH-HMM
[S2,C2,a2,gamma2,rho2] = ...
    HedgingHMM(mu_hmm,sigma_hmm,Q_hmm,T,K,r,n,put,minS,maxS,m);
display('OH-HMM done');

% OH-HMM Monte-carlo
% [S3,C3,a3,gamma3,rho3] = ...
%     HedgingHMM_MC(mu_hmm,sigma_hmm,Q_hmm,T,K,r,n,put,minS,maxS,m,N);

% OH-ARHMM
[S3,Y3,C3,a3, gamma3,rho3] = ...
    HedgingARHMM(mu_hmmvar,phi_hmmvar,sigma_hmmvar,Q_hmmvar,T,K,r,n,put,S0,m,q,N);
display('0H-ARHMM done');


%% Simulate and Hedge:

max_phi = 1;

HE_OH = 0;
CO = 0;

for i = 1:100
    
    
    R = SimHMMVAR1d(mu,phi,sigma,Q,eta0_simu,n); 
    S = S0 * exp(cumsum([0 ; R]));
    
    % Gaussian
    [HE_OH(i,1), C0(i,1), ~, ~, ~, ~,~,~,~] = ...
        Hedging_Error_Gaussian_ac(S,T,K,r,put,S1,C1,a1,rho1,max_phi);
 
    % HMM
    [HE_OH(i,2), C0(i,2), ~, ~, ~, ~, ~,~,~,~] = ...
        Hedging_Error_HMM_ac(S,T,K,r,mu_hmm,sigma_hmm,Q_hmm,eta0,put,...
            method,S2,C2,a2,rho2,max_phi);
      
    % HMM MC 
%     [HE_OH(i,3), C0(i,3), ~, ~, ~, ~, ~,~,~,~] = ...
%         Hedging_Error_HMM_ac(S,T,K,r,mu_hmm,sigma_hmm,Q_hmm,eta0,put,...
%             method,S3,C3,a3,rho3,max_phi);
        
    % ARHMM 
    [HE_OH(i,3), C0(i,3), ~, ~, ~, ~, ~,~,~,~,~] = ...
        Hedging_Error_ARHMM_ac(S,T,K,r,mu_hmmvar,phi_hmmvar,sigma_hmmvar...
            ,Q_hmmvar,eta0,put,method,S3,Y3,C3,a3,rho3,max_phi);
    
    % ARHMM MC 
%     [HE_OH(i,3), C0(i,3), ~, ~, ~, ~, ~,~,~,~] = ...
%         Hedging_Error_ARHMM_ac(S,T,K,r,mu_hmmvar,phi_hmmvar,sigma_hmmvar,
%               Q_hmmvar,eta0,put,method,S3,C3,a3,rho3,max_phi); 


        

RMSE_matrix(HE_OH)



end

nanmean(HE_OH)
RMSE_matrix(HE_OH)



%% plot density function of errors 



minX = -4
maxX =  4

bw = 0.6

figure()

[f,xi] = ksdensity(HE_OH(:,1),'width',bw);
f = f(xi>minX & xi <maxX) ; xi = xi(xi>minX & xi <maxX);
plot(xi,f,'r-x')
hold on


[f,xi] = ksdensity(HE_OH(:,2),'width',bw);
f = f(xi>minX & xi <maxX) ; xi = xi(xi>minX & xi <maxX);
plot(xi,f,'b-o')
hold on 

[f,xi] = ksdensity(HE_OH(:,3),'width',bw);
f = f(xi>minX & xi <maxX) ; xi = xi(xi>minX & xi <maxX); 
plot(xi,f,'g--o')

legend('Gaussian semi-exact','HMM semi-exact','ARHMM semi-exact','Location','NorthWest')


%%







