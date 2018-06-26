%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Compute hedging errors of ARHMM prices using optimal hedging
%
%  Input
%          
%           S: price time series
%           T: time (years) to maturity ; 1 year = 252 days ; 1 month  = 21 days 
%           K: strike price                  
%           r: annual  interest rate 
%           mu: regime means (annual values)
%           theta: autoregressive parameter - order 1
%           sigma : regime volatilities (annual values)
%           eta0: row vector of a priori probability of each regime
%           Q: transition matrix for the hidden regimes
%           put: 1 if put, 0 if call
%           m: number of points of the grid (at least 1000)
%           q: number of points on the lag return grid
%           method: 1 (most probable regime), 0 (average)
%           approx: 1 (semi-exact), 0 (monte-carlo)
%           N: Number of simulations
%           max_phi: cap on the heding position ( abs(phi) <= max_phi) )
%
%                 
%               
%      Output  
%              HE_OH: discounted hedging error ( payoff - portfolio) for
%              optimal hedging
%              phi: hedging positions
%              C0: value of the option at S0 at period 0
%              eta: (n x reg) conditional probabilities of being in regime k at 
%                   time t given observed returns r_1, ..., r_t;
%              payoff: option payoff
%              Sgrid: points of interpolation of asset price
%              Ygrid: points of interpolation of lag return
%              C: discounted value of the option at discounted prices S
%              a: discounted value of a at discounted prices S  
%              gamma: function necessary for optimal hedging (g in the
%                     paper)
%              rho: function necessary for optimal hedging (h in the paper)
% 
%
% Bruno Remillard and Massimo Caccia, April 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [HE_OH, C0, phi, eta, V, payoff, Sgrid,Ygrid,C,a,gamma,rho] = ...
            Hedging_Error_ARHMM(S,T,K,r,mu,theta,sigma,Q,eta0,put,m,q,method,approx,N,max_phi)

      if nargin == 16
          N = 10000;
      end
        
      S0 = S(1);
      Y0 = 0;
      
      n = length(S)-1;  % number of hedging periods
      reg = length(mu); % number of regimes
         
      Tp = T/n;
      rp = r*Tp;
        
      Kp = K*exp(-r*T);
        
      if approx
        [Sgrid,Ygrid,C,a,gamma,rho] = HedgingARHMM(mu,sigma,theta,Q,T,K,r,n,put,S0,m,q,N);
      else
         [Sgrid,Ygrid,C,a,gamma,rho] = HedgingARHMM_MC(mu,sigma,theta,Q,T,K,r,n,put,S0,m,q,N);
      end
   
    fprintf('Computations done for the optimal hedging strategy!\n\n');
 
    
    mu = mu*Tp ;
    sigma = sigma*sqrt(Tp);


%% Initial values


eta = zeros(n,reg);
eta(1,:) = eta0;

[X,Y] = meshgrid(Ygrid,Sgrid);

 CS = zeros(1,reg); 
 for l=1:reg  
     CS(l) = interp2D(X, Y, squeeze(C(1,:,:,l)), Y0, S0)';
 end
 
Stilde0 = S0;

R = diff(log(S));
R0 = 0 ;
    
aS = zeros(1,reg);   
phiS = zeros(1,reg); 
phi = zeros(n,1);

%%
         
        
        
         [~, reg1] = max(eta0);
         
         if method
            C0 = CS(reg1);  
         else
            C0= CS(1,:)*eta0';
         end

         V = C0;

         
          % replication starts:
          for j=1:n
              
                R1 = R(j);
                Stilde1 =  exp(-rp*j)*S(j+1);
                Delta   =  Stilde1-Stilde0;   
                
                for l=1:reg
                    aS(l) = interp2D(X,Y,squeeze(a(j,:,:,l)),R0,Stilde0);     
                    rho_ = interp1D(Ygrid,squeeze(rho(j,:,l)),R0);
                    phiS(l) = (aS(l) - V*rho_)/Stilde0 ;                              
                end
      

               if method
                   phi(j) = phiS(reg1);         
               else
                   phi(j) = phiS*eta(j,:)';
               end    

               phi(j) = max(phi(j),-max_phi); % cap on phi
               phi(j) = min(phi(j),max_phi);     
               V = V+phi*Delta  ;               
               mu_ = mu + theta * R0;    
               x = (mu_-R1).^2 ./ sigma.^2;        
               w0 = max(exp(-0.5*x ),1E-20);  
               w = w0./sigma;             
               
               z = w.*(eta(j,:)*Q)';   
               eta(j+1,:) = z/sum(z);
               [~, reg1] = max(z); 
               Stilde0 = Stilde1;     
               R0 = R1;        
              
          end
          % replication done
          
 
      
       
       if put
           payoff = max(0,Kp-Stilde1');
       else
          payoff = max(0,Stilde1'-Kp);
       end
       
       
       HE_OH =  payoff-V;
       
       

 
                   