%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Compute hedging errors of HMM prices using optimal hedging
% Grids have already been computed
%
%  Input
%          S: observed series of length n+1
%          
%           
%           T: time (years) to maturity ; 1 year = 252 days ;  
%           K: strike price                  
%           r: annual  interest rate 
%           mu: regime means (estimated from daily values)
%           sigma : regimes volatilities (daily values)
%           eta0: row vector of a priori probability of each regime
%           Q: transition matrix for the hidden regimes
%           put: 1 if put, 0 if call
%           minS: minimum value of asset price for the grid
%           maxS: maximum value of asset price for the grid
%           m: number of points of the grid
%           method: 1 (most probable regime), 0 (average)
%
%                 
%               
%      Output 
%              HE_OH: discounted hedging error ( payoff - portfolio) for
%              optimal hedging
%              HE_DH: discounted hedging error ( payoff - portfolio) for
%              the delta hedging
%              C0: value of the option at S0 at period 0
%              phi: hedging strategy 
%              eta: probabilities of each regime
%              V: discounted value of the hedging portfolio
%              Sgrid: points of interpolation of asset price
%              C: discounted value of the option at discounted prices S
%              a: discounted value of a at discounted prices S  
%              rho: function necessary for optimal hedging
     
%
% Bruno Remillard, April 5th, 2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [HE_OH, C0, phi, eta, V, payoff, Sgrid,C,a,rho] = ...
    Hedging_Error_HMM_ac(S,T,K,r,mu,sigma,Q,eta0,put,method,Sgrid,C,a,rho,max_phi)


      S0 = S(1);
      
      n = length(S)-1;   
      reg = length(mu); 
      
      %% Periodic values
       
      Tp = T/n;  
     
      rp = r*Tp;
              
      Kp = K*exp(-r*T);
   
   
    mu = mu*Tp ;
    sigma = sigma*sqrt(Tp);    


phi = zeros(n,1);

eta = zeros(n,reg);

eta(1,:) = eta0;

%% Initial values
      
 for l=1:reg
     CS(l) = interp1D(Sgrid,C(1,:,l)',S0);
 end
 
Stilde0 = S0;

R = diff(log(S));


     [~, reg1] = max(eta0);
     if method
        C0 = CS(reg1);
     else
        C0= CS(1,:)*eta0';
     end
     
     V = C0;
        
           
      for j=1:n
        
           Stilde1 =  exp(-rp*j)*S(j+1);
            Delta =  Stilde1-Stilde0;
            
            for l=1:reg
              aS(l) = interp1D(Sgrid,a(j,:,l)',Stilde0);
              phiS(l) =  (aS(l) - V*rho(j+1,l))/Stilde0;
            end
            
       if method
           phi(j) = phiS(reg1);
       else
           phi(j) = phiS*eta(j,:)';
       end

       phi(j) = max(phi(j),-max_phi);
       phi(j) = min(phi(j),max_phi);
       V = V+phi(j)*Delta;
       x = R(j);
       w0 = max(exp(-0.5*( (mu-x).^2 )./sigma.^2),1E-20);
       w = w0./sigma;
       z = w.*(eta(j,:)*Q)';
       eta(j+1,:) = z/sum(z);
       [~, reg1] = max(z);
       Stilde0 = Stilde1;
       
      end
 
      
       
       if put
           payoff = max(0,Kp-Stilde1');
       else
          payoff = max(0,Stilde1'-Kp);
       end
       
       
       HE_OH =  payoff-V;
