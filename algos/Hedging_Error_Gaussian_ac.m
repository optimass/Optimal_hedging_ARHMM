%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Compute hedging errors of geometric random walk prices using optimal hedging
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
%           N: number of simulated returns for each regime (sugg: 100000)
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
%              HE: discounted hedging error ( payoff - portfolio)
%              S: points of interpolation of asset price
%              C: discounted value of the option at discounted prices S
%              a: discounted value of a at discounted prices S  
%              rho: function necessary for optimal hedging
%    
%
% Bruno Remillard, April 5th, 2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [HE_OH, C0, phi, V, payoff, Sgrid,C,a,rho] = Hedging_Error_Gaussian_ac(S,T,K,r,put,Sgrid,C,a,rho,max_phi)

      S0 = S(1);
      
      n = length(S)-1;  % number of hedging periods
      
      
      %% Periodic values
       
      Tp = T/n;  % fraction of year by hedging period
     
      rp = r*Tp;
              
      Kp = K*exp(-r*T);


%% Initial values

phi = zeros(n,1);
      
%      C0 = interpolation_1d(S0,C(1,:)',minS,maxS);
       C0 = interp1D(Sgrid,C(1,:)',S0);

 
Stilde0 = S0;

R = diff(log(S));

     V = C0;
          
           
      for j=1:n

            Stilde1 =  exp(-rp*j)*S(j+1);
            Delta =  Stilde1-Stilde0;

              aS = interp1D(Sgrid,a(j,:)',Stilde0);
              phi(j) =  (aS - V*rho)/Stilde0;

               phi(j) = max(phi(j),-max_phi);
               phi(j) = min(phi(j),max_phi);

       V = V+phi(j)*Delta;
       Stilde0 = Stilde1;
       
       
      end

      
       
       if put
           payoff = max(0,Kp-Stilde1');
       else
          payoff = max(0,Stilde1'-Kp);
       end
       
       
       HE_OH =  payoff-V;
       

      %[HE_DH, phi_DH,V_DH,payoff_DH]=Hedging_BS(S,K,r,T,volp,put);

