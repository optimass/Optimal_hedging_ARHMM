function [S,C,a,gamma,rho] = HedgingHMM_MC(mu,vol,Q,T,K,r,n,put,minS,maxS,m,N)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% This function is an implementation by semi-exact calculations of the 
% optimal hedging solution for a put/call option, when the returns are 
% HMM Gaussian.
%
%   Input  
%       mu: average of annual returns; (reg x 1)
%       vol: volatility of annual returns; (reg x 1)
%       Q: transition matrix
%       T: time (years) to maturity ; 1 year = 252 days;1 month = 21 days; 
%       K: strike price;                  
%       r: annual interest rate; 
%       n: number of hedging periods;
%       put: 1 if put, 0 if call;
%       minS: minimum value (>0) of asset price for the grid. Note that
%             0 is always included in the computations;
%       maxS: maximum value of asset price for the grid;
%       m: number of points in the grid.
%       N: number of simulations
%
%   Output  
%       S: points of interpolation of asset price;
%       C: discounted value of the options at discounted prices S 
%          (\check C);
%       a: discounted value of a at discounted prices S  (\check a);
%       gamma: function necessary for optimal hedging
%       rho  : function necessary for optimal hedging
%
%
% 
% Bruno Remillard and Massimo Caccia, April 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


       
      %% check-up
      
      reg = length(mu);
      if ~( length(vol)*length(Q) ==  reg^2 )
          error('problem with number of regimes')
      end    
       
      %% periodic values 
      Tp = T/n;
      rp = r*Tp;
   
      Kp = K*exp(-r*T);
      dt = (maxS-minS)/(m-1); 
   
      %%
      S1  = (minS:dt:maxS)';
      S = [0;S1];
      
      volp = vol*sqrt(Tp);  


     mu2 = mu*Tp -rp;  
     sigma2 = volp.^2; 
    
     %% MC sample

     % old way :
%      X = repmat(mu2',N,1) + randn(N,reg) .* repmat(volp',N,1); 
     
     % new way
     X = randn(N,reg);
     X = ( X - mean(X)) ./ std(X);
     X = X .* repmat(volp',N,1) + repmat(mu2',N,1) ;
     
     M1 = mean(exp(X)-1);
     M2 = mean((exp(X)-1).^2);
     
     

   
  %% 
  

  
  gamma = zeros(n+1,reg);
  rho   = zeros(n+1,reg);
  a     = zeros(n,reg); 
  b     = zeros(n,reg);
  
  gamma(n+1,:) = 1;
  
  for k = n:-1:1
      
      if k == 9
          1;
      end
      
     b(k,:)     = Q * (gamma(k+1,:)' .* M1');
     a(k,:)     = Q * (gamma(k+1,:)' .* M2');
     rho(k+1,:) = b(k,:) ./ a(k,:);
     gamma(k,:) = Q*gamma(k+1,:)' - (rho(k+1,:) .* b(k,:))' ;
     
  end
  
g = a;
  
  
  %% PRE PRICING
      
   C = zeros((n+1),(m+1),reg);
   a = zeros(n,(m+1),reg); 
             

   if put
       C(n+1,:,:)  = repmat( max(Kp-S,0) , 1 , reg );    % put
   else 
       C(n+1,:,:)  = repmat( max(S-Kp,0) , 1 , reg ) ;   % call
   end



  C(1:n,1,:) =  squeeze(C(n+1,1,:))' .* gamma(1:end-1,:) ; 
  a(:,1,:) = reshape(C(1:n,1,:),[n,reg]) .* rho(2:end,:);

  C_next = zeros(N,reg);

  
   %% PRICING STARTS 
   
for k = n:-1:1
   
    
    fprintf(' Period  %.2f \n', k-1);
      
           
    
    %%   C and a 
        
         for p = 2:(m+1)
           
           s=S(p); 
           
           s_next = s*exp(X);
           
           for j = 1:reg
              C_next(:,j) = interp1(S,C(k+1,:,j),s_next(:,j),'linear','extrap');
           end
             
           for l = 1:reg 
              
                 
               % compute C:  
               C(k,p,l) = mean( ( 1 - rho(k+1,l).*(exp(X)-1) ) .* C_next  )* Q(l,:)' ;

               % compute a:
               temp_a = mean( (exp(X)-1) .* C_next ) * Q(l,:)' ;
               a(k,p,l) = temp_a  / g(k,l) ;


            
            end % end of regime loop
         end % end of price loop        
end % end of period loop


    temp_gamma = repmat(gamma,1,1,m+1);
    temp_gamma = permute(temp_gamma,[1 3 2]);
    C = C ./ temp_gamma;

    a(:,1,:) = squeeze(a(:,1,:)) ./ gamma(1:end-1,:);
    
    
   
