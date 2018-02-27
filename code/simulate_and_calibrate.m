%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% In this script, we will simulate and calibrate the following models:
%
% 1) HMM: univariate Gaussian hidden Markov model
%       EstHMM for calibration
%       GofHMM for Goodness-of-fit test
%
% 2) ARHMM: univariate autoregressive Gaussian hidden Markov model
%       EstARHMM for calibration
%       GofARHMM for Goodness-of-fit test
%
% 3) VHMM: multivariate Gaussian hidden Markov model
%       EstVHMM for calibration
%       GofVARHMM for Goodness-of-fit test
%
% 4) VARHMM: multivariate autoregressive Gaussian hidden Markov model
%       EstVHMM for calibration
%       GofVARHMM for Goodness-of-fit test
%
% models from the paper :
%       "Option Pricing and Hedging for Discrete 
%        Time Autoregressive Hidden Markov Model"
%
% https://papers.ssrn.com/sol3/papers.cfm?abstract_id=2995944
%
%
%%%%%%%%%%%%%%%%

%% 0 :: Set global parameters

reg = 3;
d=2;
max_iter = 1000;
prec = 0.0001;
N = 1000;
n = 1*252;
T = n/252;
r = 0.01;
eta0 = 1;


%% 1 :: HMM

% set params
mu = [0.01 0.02 -0.03];
sigma = [0.02 0.04 0.08];
Q = [0.95  0.025 0.025 ; 
     0.025 0.95  0.025 ;
     0.025 0.025 0.95];
 
% simulate
[y,~] = SimHMM(mu,sigma,Q,eta0,n);

% calibrate
[mu_,sigma_,Q_,eta,nu,Z] = EstHMM(y,reg,max_iter,eps);

% goodness-of-fit
out_hmm = GofHMM(y,reg,max_iter,prec,N)

%% 2 :: ARHMM

% set params
mu = [0.01 0.02 -0.03];
phi = [0.2 0 -0.2];
sigma = [0.02 0.04 0.08];
Q = [0.95  0.025 0.025 ; 
     0.025 0.95  0.025 ;
     0.025 0.025 0.95];
 
% simulate
[y,~] = SimARHMM(mu,phi,sigma,Q,eta0,n);

% calibrate
[mu_,phi_,sigma_,Q_,eta,nu,Z] = EstARHMM(y,reg,max_iter,eps);

% goodness-of-fit
out_arhmm = GofARHMM(y,reg,max_iter,prec,N)


%% 3 :: VHMM
% might need a bigger time horizon for more accurate results
% e.g. 
n = 10*252;
T = n/252;
N = 100;

% set params
mu = [0.01 0.02  -0.03;
      0.05 -0.04 0.00];
A = zeros(d,d,reg); % note that A is covariance (not volatility)
A(:,:,1) = [0.05 0.03 ; 0.03 0.03];
A(:,:,2) = [0.1 -0.05 ; -0.05 0.12];
A(:,:,3) = [0.05 0.00 ; 0.00 0.02];

Q = [0.95  0.025 0.025 ; 
     0.025 0.95  0.025 ;
     0.025 0.025 0.95];
 
% simulate
[y,~] = SimVHMM(mu,A,Q,eta0,n);
 
% calibrate
[mu_,A_,Q_,eta,nu,Z] = EstVHMM(y,reg,max_iter,eps);
 
% goodness-of-fit
out_vhmm = GofVHMM(y,reg,max_iter,prec,N);


%% 4 :: VARHMM

% set params
mu = [0.01 0.02  -0.03;
      0.05 -0.04 0.00];
phi = zeros(d,d,reg);
phi(:,:,1) = [0.0  0.3  ; 0.0 0.0];
phi(:,:,2) = [0.5  0.0  ; 0.0 0.0];
phi(:,:,3) = [0.0  0.0 ; -0.4 0.0];
A = zeros(d,d,reg);
A(:,:,1) = [0.05 0.03 ; 0.03 0.03];
A(:,:,2) = [0.1 -0.05 ; -0.05 0.12];
A(:,:,3) = [0.05 0.00 ; 0.00 0.02];
Q = [0.95  0.025 0.025 ; 
     0.025 0.95  0.025 ;
     0.025 0.025 0.95];
 
 % simulate
[y,~] = SimVARHMM(mu,phi,A,Q,eta0,n);
 
% calibrate
[mu_,phi_,A_,Q_,eta,nu,Z] = EstVARHMM(y,reg,max_iter,eps);
 
% goodness-of-fit
out_varhmm = GofVARHMM(y,reg,max_iter,prec,N);
 

















