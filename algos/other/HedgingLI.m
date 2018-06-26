function [S,C,a,c,rho,phi1] = HedgingLI(mu,vol,T,K,r,n,put,minS,maxS,S0,m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% This function is an implementation by semi-exact calculations of the 
% optimal hedging solution for a put/call option, when the returns are 
% Gaussian. The notations are taken from Example 3.4.7. 
%
%   Input  
%       mu: average of annual returns;
%       sigma: volatility of annual returns;
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
%       c: function necessary for optimal hedging
%       rho  : function necessary for optimal hedging
%       phi1 : initial allocation of the risky asset 
%
%
% Example: [S,C,a,c,rho,phi1] = HedgingLI(0.09,0.06,1,100,0.05,22,0,80,120,100,2000);
% 
% Bruno Remillard, April 5th, 2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
     
       
      %% periodic values 
      Tp = T/n;
      rp = r*Tp;
   
      Kp = K*exp(-r*T);
      dt = (maxS-minS)/(m-1); 
   
      %%
      S1  = (minS:dt:maxS)';
      S = [0;S1];
      
      volp = vol*sqrt(Tp);  % periodic vol
 
     I0 =  zeros((m+1),(m+1));  
     I1 =  I0;
     I2 =  I0;
 
     mu2 = mu*Tp -rp;  %periodic mean of excess returns
     sigma2 = volp^2;
     c1 = exp(mu2+0.5*sigma2);
     c2 = exp(2.0*(mu2+sigma2));
     M1 = c1-1.0;
     M2 = c2+1.0 -2.0*c1;
       
     I0(1,1) = 1;
     I1(1,1) = c1;
     I2(1,1) = c2;

     ind2 =(2:(m+1))';
       
 %%Equations for A and B at q=0 and q=m 

         s=S1(1);
       
           d0 = (log(s./S1)-mu2)/volp;
           d1 = d0 - volp;
           d2 = d1 - volp;
           I0(1,ind2) = normcdf(d0);
           I1(1,ind2) = normcdf(d1)*c1;
           I2(1,ind2) = normcdf(d2)*c2;

           s= S1(m);
           d0 = (log(s./S1)-mu2)/volp;
           d1 = d0 - volp;
           d2 = d1 - volp;
           I0(m+1,ind2) = 1.0-normcdf(d0);
           I1(m+1,ind2) = (1.0-normcdf(d1))*c1;
           I2(m+1,ind2) = (1.0-normcdf(d2))*c2;

   
   for k=2:m
          
           s = S(k+1);
           d0 = (log(s./S1)-mu2)/volp;
           d1 = d0 - volp;
           d2 = d1 - volp;
           I0(k,ind2) = normcdf(d0);
           I1(k,ind2) = normcdf(d1)*c1;
           I2(k,ind2) = normcdf(d2)*c2;

           s = S(k);
           d0 = (log(s./S1)-mu2)/volp;
           d1 = d0 - volp;
           d2 = d1 - volp;
           I0(k,ind2) =  I0(k,ind2)-normcdf(d0');
           I1(k,ind2) =  I1(k,ind2)-normcdf(d1')*c1;
           I2(k,ind2) =  I2(k,ind2)-normcdf(d2')*c2;
         
   end
        
      
   rho = M1/M2; 
   c = (1-rho*M1);
   C = zeros((n+1),(m+1));
   a = zeros(n,(m+1));
             

   C(n+1,:)  = max(S-Kp',0);   % call
   if(put)
       C(n+1,:)  = max(Kp-S',0);
   end


  C(1:n,1) = C(n+1,1);
  a(:,1) = C(1:n,1)*rho;
  
   %%
 A = zeros(n,(m+1)) ; 
  B =   A;

for k = (n:-1:1)
   
    
    fprintf(' Period  %d \n', k-1);
      
       
            
    %%   A and B 
     A(k,1:m) = C(k+1,1:m); 
     B(k,1:m) = (C(k+1,2:(m+1))-A(k,1:m))./(S(2:m+1)-S(1:m))';   
                   
     A(k,1:m) = A(k,1:m)-B(k,1:m).*S(1:m)'; 
       
     A(k,m+1) = A(k,m); 
     B(k,m+1) = B(k,m);      
    
%%    C and a
        
         for i=2:(m+1)
            
           s=S(i);
           sumC = A(k,:)*( (1.0+rho)*I0(:,i)-rho*I1(:,i)) +  s*B(k,:)*( (1.0+rho)*I1(:,i) - rho*I2(:,i));
           suma = A(k,:)*(I1(:,i)-I0(:,i)) +s*B(k,:)*(I2(:,i) -I1(:,i));
                 
                      
            a(k,i) = suma/M2;
            C(k,i) = sumC/c;
         
         end
                 
     

end
     AA(1:m) = C(1,1:m); 
     BB(1:m) = (C(1,2:(m+1))-AA(1:m))./(S(2:m+1)-S(1:m))';   
                     
                  
     AA(1:m) = AA(1:m)-BB(1:m).*S(1:m)'; 
            
     AA(m+1) = AA(m);
     BB(m+1) =BB(m);  
     
     %% phi1
       phi1 = zeros(m,1);
    
       z = (a(1,2:end) - C(1,2:end)*rho)';
       phi1(:) = z./S1;
       
        
     
     
%        CS0 = interpolation_new_1d(S0,C(1,:)',minS,maxS);
%        aS0 = interpolation_new_1d(S0,a(1,:)',minS,maxS);
       
      
       CS0 = interp1(S,C(1,:)',S0);
       aS0 = interp1(S,a(1,:)',S0);
       
       phiS0 =  (aS0 - CS0*rho)/S0;
  
   [CS0 phiS0]
     
    addpath('../Chapter01')
    [C_BS, P_BS] = FormulaBS(S0,K,r,T,vol)
    [deltaC_BS,deltaP_BS] = FormulaBSGreeks(S0,K,r,T,vol)

     
    [C_DH, C1_DH] = FormulaBS(S1,K,r,T,vol);
    [phi_DH,phi1_DH] = FormulaBSGreeks(S1,K,r,T,vol);
    rmpath('../Chapter01')
      if(put)
       phi_DH = phi1_DH;
         C_DH = C1_DH;
      end

                      
       figure
       plot(S1,C(1,2:end),'r--')
       hold on
       plot(S1,C_DH)
       if(put)
          title('Put values for Black-Sholes model')
        else
           title('Call values for Black-Sholes model')
         
       end
       legend('Delta hedging','Optimal hedging','Location','NorthWest');

       figure
       plot(S1,phi1,'r--')
       hold on
       plot(S1,phi_DH)
       if(put)
          title('Initial investment in the put for Black-Sholes model')
         else
          title('Initial investment in the call for Black-Sholes model')
        end
       legend('Delta hedging','Optimal hedging','Location','NorthWest');

     
       
