function [S,Y,C,a,g,h] = HedgingARHMM_MC(mu,phi,vol,Q,T,K,r,n,put,S0,m,q,N)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% This function is an implementation by semi-exact calculations of the 
% optimal hedging solution for a put/call option, when the returns are 
% HMM Gaussian.
%
%   Input  
%       mu: average of annual returns; (reg x 1)
%       phi: autoregressive parameter - order 1 (reg x 1)
%       vol: volatility of annual returns; (reg x 1)
%       Q: transition matrix
%       T: time (years) to maturity ; 1 year = 252 days;1 month = 21 days; 
%       K: strike price;                  
%       r: annual interest rate; 
%       n: number of hedging periods;
%       put: 1 if put, 0 if call;
%       S0: starting price
%       m: number of points in the asset price grid.
%       q: number of points on the lag return grid
%       N: number of simulations
%
%   Output  
%       S: points of interpolation of asset price;
%       C: discounted value of the options at discounted prices S
%       a: discounted value of a at discounted prices S 
%       gamma (g): function necessary for optimal hedging
%       rho (h): function necessary for optimal hedging
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
      
      if q<20
          error('q has to be bigger  >= 20')
      end  
       
      %% periodic values 
      Tp = T/n;
      rp = r*Tp;
   
      Kp = K*exp(-r*T); 
      
      volp = vol*sqrt(Tp); 
      mup = mu*Tp ;       
      mu2 = mu*Tp -rp' ; 
   
      %%
      
      %%% trick for S grid %%%%
      for i = 1:N
        R(i,:) = SimHMMVAR1d(mu2,phi,volp,Q,1,n); 
      end
      R = S0 .* exp(cumsum(R,2));
      R = R(:);
      minS = min(R)-0.005 ; maxS = max(R)+0.005 ; 
      S = [minS; prctile(R,0.01) ; prctile(R,0.1) ; prctile(R,1) ; prctile(R,2.5) ] ;
      for i = 1:m-10
          S = [S ; prctile(R,i/(m-9)*100)];
      end
      S = [S ; prctile(R,97.5) ; prctile(R,99) ; prctile(R,99.9) ; prctile(R,99.99); maxS];
      S = sort(S);
      S = [0 ; S];
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      %%% trick for Y grid %%%%
      R = SimHMMVAR1d(mup,phi,volp,Q,1,100000); 
      minY = min(R)-0.005 ; maxY = max(R)+0.005 ; 
      Y = [minY; prctile(R,0.01) ; prctile(R,0.1) ; prctile(R,1) ; prctile(R,2.5) ] ;
      for i = 1:q-10
          Y = [Y ; prctile(R,i/(q-9)*100)];
      end
      Y = [Y ; prctile(R,97.5) ; prctile(R,99) ; prctile(R,99.9) ; prctile(R,99.99); maxY];
      Y = sort(Y);
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
     %% MC sample
     

% Old Way :
%      X = randn(N,reg);
%      X = ( X - repmat(mean(X),N,1) ) ./ repmat(std(X),N,1);
%      X = X .* repmat(volp',N,1) + repmat( (1-phi').*mup' ,N,1) ; %%%%%
%      X = repmat(X,1,1,q);   
%      lag = phi*Y';
%      lag = repmat(lag,1,1,N);
%      lag = permute(lag,[3 1 2]);
%      X = X + lag;
%      X2 = X - rp;
     
     X = randn(N,reg);
     X = ( X - repmat(mean(X),N,1) ) ./ repmat(std(X),N,1);
     
     X = X .* repmat(volp',N,1) + repmat( (1-phi').*mup' ,N,1) ; %%%%%
     X = repmat(X,1,1,q);   
     lag = phi*Y';
     lag = repmat(lag,1,1,N);
     lag = permute(lag,[3 1 2]);
     X = X + lag;
     X2 = X - rp;  
    
   
  %% 
  
  % memoire notation
  
  g = zeros(n+1,q,reg);
  a = zeros(n,q,reg); 
  b = zeros(n,q,reg);
  h = zeros(n,q,reg);
  
  g(n+1,:,:) = 1;
  
  g_next = zeros(N,reg);

  
  for k = n:-1:1
       
    for v = 1:q
        
        for i = 1:reg
            g_next(:,i) = interp1(Y,g(k+1,:,i),X(:,i,v),'linear','extrap');
        end
          
        
        for i = 1:reg
            
            a(k,v,i) = mean(( (exp(X2(:,:,v))-1).^2 .* g_next)) * Q(i,:)' ;
            b(k,v,i) = mean(( (exp(X2(:,:,v))-1)    .* g_next)) * Q(i,:)' ;
                               
            h(k,v,i) = b(k,v,i) / a(k,v,i) ;
            g(k,v,i) = mean(g_next * Q(i,:)') - h(k,v,i)*b(k,v,i);
        
        end
    end
    
  end
   
  %% PRE PRICING
      
   C = zeros((n+1),(m+1),q,reg);
   A = zeros(n,(m+1),q,reg); 
             

   if put
       C(n+1,:,:,:)  = repmat( max(Kp-S,0) , 1 , q, reg );    % put
   else 
       C(n+1,:,:,:)  = repmat( max(S-Kp,0) , 1 , q, reg ) ;   % call
   end

   
   
  % book page 79 and BrunoPdf --> fill when s=0
  temp = repmat(squeeze(C(n+1,1,:,:)),1,1,n);
  temp = permute(temp,[3 1 2]);
  C(1:n,1,:,:) = temp ;  

  A(:,1,:,:) = squeeze(C(1:n,1,:,:)) .* h;

  C_next = zeros(N,reg);

   %% PRICING STARTS 
   
for k = n:-1:1
   
    
    fprintf(' Period  %.2f \n', k-1);
      
           
    
    %%   C and a 
        
         for p = 2:(m+1)
           
           s=S(p); 
           
           for u = 1:q
           
             
               s_next = s*exp(X2(:,:,u));

               for i = 1:reg
                  C_next(:,i) = interp1(S, C(k+1,:,u,i), s_next(:,i), 'linear', 'extrap');
                  g_next(:,i) = interp1(Y, g(k+1,:,i)  , X(:,i,u)   , 'linear', 'extrap');
               end

               for i = 1:reg 
                    
                   % compute C:              
                    C(k,p,u,i) = mean( ( 1 - h(k,u,i).*(exp(X2(:,:,u))-1) ) .* C_next .* g_next) * Q(i,:)' ;
                    C(k,p,u,i) = C(k,p,u,i) / g(k,u,i);
                    
                   % compute a: 
                    A(k,p,u,i) = mean( (exp(X2(:,:,u))-1) .* C_next .* g_next ) * Q(i,:)' ; 

            
                end % end of regime loop
            end % end of lag return loop
         end % end of price loop       
end % end of period loop

     
     temp_a = repmat(a,1,1,1,m);
     temp_a = permute(temp_a,[1 4 2 3]);
     a= A;
     a(:,2:end,:,:) = a(:,2:end,:,:) ./ temp_a;

     
