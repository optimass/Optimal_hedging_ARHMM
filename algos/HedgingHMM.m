function [S,C,a,gamma,rho] = HedgingHMM(mu,vol,Q,T,K,r,n,put,minS,maxS,m)
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
% Bruno Remillard and Massimo Caccia, April 15, 2016
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
     mup = mu*Tp -rp;  
     sigma2 = volp.^2; 
     c0 = 1;
     c1 = exp(mu2+0.5*sigma2);  
     c2 = exp(2.0*(mu2+sigma2));  
     M1 = c1-1.0;   
     M2 = c2+1.0 -2.0*c1;  
       
      %% create grids
      
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       c = c0;
       theta = 0;
       I0 = create_grid3D(m,reg,S1,c,theta,mup,volp);           
       c = c1;
       theta = 1;
       I1 = create_grid3D(m,reg,S1,c,theta,mup,volp);
       c = c2;
       theta = 2;
       I2 = create_grid3D(m,reg,S1,c,theta,mup,volp);
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   
  %% 
  
  
  gamma = zeros(n+1,reg);
  rho   = zeros(n+1,reg);
  a     = zeros(n,reg); 
  b     = zeros(n,reg);
  
  gamma(n+1,:) = 1;
  
  for k = n:-1:1
     b(k,:)     = Q * (gamma(k+1,:)' .* M1);
     a(k,:)     = Q * (gamma(k+1,:)' .* M2);
     rho(k+1,:) = b(k,:) ./ a(k,:);
     gamma(k,:) = Q*gamma(k+1,:)' - (rho(k+1,:) .* b(k,:))';
     
  end
  

g = a;

  %% PRE PRICING
      
   C = zeros((n+1),(m+1),reg);
   a = zeros(n,(m+1),reg); 
             

       C(n+1,:,:)  = repmat( max(S-Kp,0) , 1 , reg ) ;   % call
   if put
       C(n+1,:,:)  = repmat( max(Kp-S,0) , 1 , reg );    % put
   end


  C(1:n,1,:) =  squeeze(C(n+1,1,:))' .* gamma(1:end-1,:) ;  
  a(:,1,:) = squeeze(C(1:n,1,:)) .* rho(2:end,:); 
  
  A = zeros(n,(m+1),reg) ; 
  B = A;

  
   %% PRICING STARTS 
   
for k = n:-1:1
   
    
    fprintf(' Period  %.2f \n', k-1);
      
       
            
    %%   A and B 
  
    
     ds = repmat((S(2:m+1)-S(1:m))',1,1,reg); 
     B(k,1:m,:) = ( C(k+1,2:(m+1),:)-C(k+1,1:m,:) ) ./ ds;
     
     temp = repmat(S(1:m)',1,1,reg);                 
     A(k,1:m,:) = C(k+1,1:m,:)-B(k,1:m,:).*temp; 
       
     A(k,m+1,:) = A(k,m,:); 
     B(k,m+1,:) = B(k,m,:);
    
    %%   C and a 
        
         for p = 2:(m+1)
           
           s=S(p);  
             
            for l = 1:reg 
           
           
           % first part of the sum:
           Csum0 = sum(squeeze(A(k,:,:)) .* squeeze(I0(:,p,:)))'  ;         
           Csum0 = (1+rho(k+1,l)) * (Q(l,:) * Csum0) ;

           asum0 = sum(squeeze(A(k,:,:)) .* squeeze(I0(:,p,:)))';
           asum0 = -1 * (Q(l,:) * asum0);
           
           % second part of the sum 
           
           Csum1 = s * (1+rho(k+1,l)) * squeeze(B(k,:,:));          
           Csum1 = Csum1 - rho(k+1,l) * squeeze(A(k,:,:));        
           Csum1 = Csum1 .* squeeze(I1(:,p,:));                 
           Csum1 = sum(Csum1)';                   
           Csum1 = Q(l,:) * Csum1;                             
           
           asum1 = squeeze(A(k,:,:)) - s * squeeze(B(k,:,:));      
           asum1 = asum1 .* squeeze(I1(:,p,:));
           asum1 = sum(asum1)';
           asum1 = Q(l,:)* asum1;
           
           % third part of the sum
           
           Csum2 = squeeze(B(k,:,:)) .* squeeze(I2(:,p,:));          
           Csum2 = sum(Csum2)';                                    
           Csum2 = Q(l,:) * Csum2;                                  
           Csum2 = -1 * s * rho(k+1,l) * Csum2  ;                  
           
           asum2 = squeeze(B(k,:,:)) .* squeeze(I2(:,p,:));      
           asum2 = sum(asum2)'; 
           asum2 = Q(l,:) * asum2;  
           asum2 = s * asum2;
                      
           % sum everything up:
           
           sumC     = Csum0 + Csum1 + Csum2;
           C(k,p,l) = sumC ;
           
           suma     = asum0 + asum1 +asum2;
           a(k,p,l) = suma / g(k,l) ;

            
            end % end of regime loop
         end % end of price loop         
end % end of period loop


    temp_gamma = repmat(gamma,1,1,m+1);
    temp_gamma = permute(temp_gamma,[1 3 2]);
    C = C ./ temp_gamma;

    a(:,1,:) = squeeze(a(:,1,:)) ./ gamma(1:end-1,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
     
     
     
     
function  [grid] = create_grid3D(m,reg,S1,c,theta,mup,volp)

grid = nan(m+1,m+1,reg);

grid(1,1,:) = c;  % its the expectation on the lowest boarder, where s = 0

% useful for computation
mup_  = repmat(mup',m,1);
volp_ = repmat(volp',m,1);
ind2 =(2:(m+1))'; % all indice except first
c_   = repmat(c',m,1);

%% fill top and bottom borders (at q=0 and q=m)

s  = S1(1); % prix le plus bas (sauf 0 )
x  = log(s./S1);
x  = repmat(x,1,reg);
d = (x-mup_)./ volp_ - theta*volp_;
grid(1,ind2,:) = normcdf(d).*c_;

s  = S1(m); % prix le plus élevé
x  = log(s./S1);
x  = repmat(x,1,reg);
d = (x-mup_)./ volp_ - theta*volp_;
uno = ones(length(ind2),reg);
grid(m+1,ind2,:) = (uno-normcdf(d)).*c_;

 %% on calcule les autre
 
for k=2:m
          
           % little trick to calculate the interval :
           s  = S1(k);
           x  = log(s./S1);
           x  = repmat(x,1,reg);
           d = (x-mup_)./ volp_ - theta*volp_;
           grid(k,ind2,:) = normcdf(d).*c_;
           
           grid0 = grid;
           
           s = S1(k-1);
           x = log(s./S1);
           x  = repmat(x,1,reg);
           d = (x-mup_)./ volp_ - theta*volp_;  
           grid(k,ind2,:) = squeeze(grid0(k,ind2,:)) - normcdf(d) .* c_;

end



         
         
     

