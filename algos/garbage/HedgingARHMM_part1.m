function [out] = HedgingARHMM(mu,phi,vol,Q,T,K,r,n,put,S0,m,q,N)
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
%
%   Output  
%       S: points of interpolation of asset price;
%       Y: points of interpolation of lag return;
%       C: discounted value of the options at discounted prices S
%       a: discounted value of a at discounted prices S 
%       g: function necessary for optimal hedging
%       h: function necessary for optimal hedging
%
%
% 
% Bruno Remillard and Massimo Caccia, April 15, 2017
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
      
      volp = vol*sqrt(Tp);  
      mup = mu*Tp ;      
%       mu2 = mu*Tp -rp' ;  
      mu2 = mu*Tp ;
      sigma2 = volp.^2; 
   

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

      
      %% create grids for Y
      
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       theta = 0;
       Y_M0 = create_grid3D_M(q,reg,Y,theta,mu2,phi,volp);           
       theta = 1;
       Y_M1 = create_grid3D_M(q,reg,Y,theta,mu2,phi,volp);
       theta = 2;
       Y_M2 = create_grid3D_M(q,reg,Y,theta,mu2,phi,volp);
       %%%%%%%%%%%%%%%%%%
       theta = 0;
       Y_Mp0 = create_grid3D_Mprime(q,reg,Y,theta,mu2,phi,volp,Y_M0);           
       theta = 1;
       Y_Mp1 = create_grid3D_Mprime(q,reg,Y,theta,mu2,phi,volp,Y_M1);
       theta = 2;
       Y_Mp2 = create_grid3D_Mprime(q,reg,Y,theta,mu2,phi,volp,Y_M2); 
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% in order to not squeeze

    Y_M0 = permute(Y_M0,[1 3 2]);
    Y_M1 = permute(Y_M1,[1 3 2]);
    Y_M2 = permute(Y_M2,[1 3 2]);
    
    Y_Mp0 = permute(Y_Mp0,[1 3 2]);
    Y_Mp1 = permute(Y_Mp1,[1 3 2]);
    Y_Mp2 = permute(Y_Mp2,[1 3 2]);
       
       
  %% 
  
  % memoire notation
  
  g = zeros(n+1,q,reg);
  a = zeros(n,q,reg); 
  b = zeros(n,q,reg);
  h = zeros(n,q,reg);
  
  g(n+1,:,:) = 1;
 
  
  for k = n:-1:1
      
      F = squeeze(g(k+1,:,:));
      [A_g,B_g] = create_A_B_2D(Y,F,q,reg) ;
         
        for v = 1:q
                for i = 1:reg
            
                    b(k,v,i) = sum(- (A_g.*Y_M0(:,:,v) + B_g.*Y_Mp0(:,:,v)) ...
                                   + exp(-rp)*A_g.*Y_M1(:,:,v) + exp(-rp)*B_g.*Y_Mp1(:,:,v) ) * Q(i,:)' ;

                    a(k,v,i) = sum(A_g.*Y_M0(:,:,v) + B_g.*Y_Mp0(:,:,v) ...
                                    - 2*exp(-rp)*A_g.*Y_M1(:,:,v) - 2*exp(-rp)*B_g.*Y_Mp1(:,:,v) ...
                                        + exp(-2*rp)*A_g.*Y_M2(:,:,v) + exp(-2*rp)*B_g.*Y_Mp2(:,:,v) ) * Q(i,:)' ;

                    h(k,v,i) = b(k,v,i) / a(k,v,i) ;
                    g(k,v,i) = sum( A_g.*Y_M0(:,:,v) + B_g.*Y_Mp0(:,:,v) ) * Q(i,:)' - h(k,v,i)*b(k,v,i);
                    
                end
        end
  end

  %% create grids for C 
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  theta = 0;
  SY_M0 = create_grid5D_M(m,q,reg,S,Y,theta,mu2,phi,volp);
  theta = 1;
  SY_M1 = create_grid5D_M(m,q,reg,S,Y,theta,mu2,phi,volp);
  theta = 2;
  SY_M2 = create_grid5D_M(m,q,reg,S,Y,theta,mu2,phi,volp);
  %%%%%%%%%%%%%%%%%%%%%
  theta = 0;
  SY_Mp0 = create_grid5D_Mprime(m,q,reg,S,Y,theta,mu2,phi,volp,SY_M0);
  theta = 1;
  SY_Mp1 = create_grid5D_Mprime(m,q,reg,S,Y,theta,mu2,phi,volp,SY_M1);
  theta = 2;
  SY_Mp2 = create_grid5D_Mprime(m,q,reg,S,Y,theta,mu2,phi,volp,SY_M2);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  
  %% in order to not squeeze

    SY_M0 = permute(SY_M0,[3 4 5 1 2]);
    SY_M1 = permute(SY_M1,[3 4 5 1 2]);
    SY_M2 = permute(SY_M2,[3 4 5 1 2]);
    
    SY_Mp0 = permute(SY_Mp0,[3 4 5 1 2]);
    SY_Mp1 = permute(SY_Mp1,[3 4 5 1 2]);
    SY_Mp2 = permute(SY_Mp2,[3 4 5 1 2]);
    
  %%
  
  out.S = S;
  out.Y = Y;
  
  out.a = a;
  out.b = b;
  out.g = g;
  out.h = h;
  
  out.SY_M0 = SY_M0;
  out.SY_M1 = SY_M1;
  out.SY_M2 = SY_M2;
  out.SY_Mp0 = SY_Mp0;
  out.SY_Mp1 = SY_Mp1;
  out.SY_Mp2 = SY_Mp2;
  
  

   
  %% PRE PRICING
  
  
%    % \check C case:    
%    C = zeros((n+1),(m+1),q,reg);
%    A = zeros(n,(m+1),q,reg); 
%              
% 
%    if put
%        C(n+1,:,:,:)  = repmat( max(Kp-S,0) , 1 , q, reg );    % put
%    else 
%        C(n+1,:,:,:)  = repmat( max(S-Kp,0) , 1 , q, reg ) ;   % call
%    end
% 
% 
%    
%   % book page 79 and BrunoPdf --> fill when s=0
%   temp = repmat(squeeze(C(n+1,1,:,:)),1,1,n);
%   temp = permute(temp,[3 1 2]);
%   C(1:n,1,:,:) = temp .* g(1:end-1,:,:);   
% 
%   A(:,1,:,:) = squeeze(C(1:n,1,:,:)) .* h;


% 
%    %% PRICING STARTS 
%    
% for k = n:-1:1
%    
%     fprintf(' Period  %.2f \n', k-1);
%       
%     for p = 2:(m+1)
%            
%         s_q=S(p); 
%            
%             for u = 1:q
%                 
%                 F = squeeze(C(k+1,:,:,:));
%                [A_C,B1_C,B2_C,B3_C] = create_A_B123_3D(S,Y,F,m,q,reg) ; 
%                
%                
%                for i = 1:reg 
%                    
% C(k,p,u,i) =  Q(i,:) * squeeze(sum(sum( (1+h(k,u,i)) .* A_C .* SY_M0(:,:,:,p,u) ...  
%                     + (1+h(k,u,i)) .* B2_C .* SY_Mp0(:,:,:,p,u) ...
%                     + ( s_q*(1+h(k,u,i)) .* B1_C - h(k,u,i) .* A_C ) .*  SY_M1(:,:,:,p,u) * exp(-rp) ... 
%                     + ( s_q*(1+h(k,u,i)) .* B3_C - h(k,u,i) .* B2_C) .* SY_Mp1(:,:,:,p,u) * exp(-rp)...
%                     - h(k,u,i)*s_q .* B1_C .*  SY_M2(:,:,:,p,u) * exp(-2*rp)...                      
%                     - h(k,u,i)*s_q .* B3_C .* SY_Mp2(:,:,:,p,u) * exp(-2*rp)))) ;
% 
% A(k,p,u,i) = Q(i,:) * squeeze(sum(sum( - A_C .* SY_M0(:,:,:,p,u) ...    
%                     - B2_C .* SY_Mp0(:,:,:,p,u) ...
%                     + ( A_C - s_q.*B1_C  ) .* SY_M1(:,:,:,p,u) * exp(-rp)...
%                     + ( B2_C - s_q.*B3_C ) .* SY_Mp1(:,:,:,p,u) * exp(-rp) ...
%                     + s_q.*B1_C .* SY_M2(:,:,:,p,u) * exp(-2*rp)  ...  
%                     + s_q.*B3_C .* SY_Mp2(:,:,:,p,u) * exp(-2*rp) ))) ;
% 
%             
%                 end % end of regime loop
%                 
%             end % end of lag return loop
%             
%     end % end of price loop
%                  
% end % end of period loop
% 
% 
%     temp_g = repmat(g,1,1,1,m+1);
%     temp_g = permute(temp_g,[1 4 2 3]);
%     C = C ./ temp_g;
% 
%     A(:,1,:,:) = squeeze(A(:,1,:,:)) ./ g(1:end-1,:,:);
%      
% 
%      temp_a = repmat(a,1,1,1,m);
%      temp_a = permute(temp_a,[1 4 2 3]);
%      a = A;
%      a(:,2:end,:,:) = a(:,2:end,:,:) ./ temp_a;  
%      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
     
     
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     function  [grid] = create_grid3D_M(m,reg,Y,theta,mup,phi,volp)

grid = zeros(m+1,m,reg);

%% add lag return dim to mu

       mup_ = (1-phi').*mup';
       lag = Y*phi';
       mup_ = lag+mup_;
       
       volp_ = repmat(volp',m,1);
       
      
%% compute constant c (in front of the intervals)

      c = exp(theta*mup_ + theta^2 * volp_.^2 /2);
       
%% fill top and bottom borders (at q=0 and q=m)

    x  = Y(1);
    d = (x-mup_)./ volp_ - theta*volp_; % Kappa in memoire
    grid(1,:,:) = normcdf(d).*c;

    x  = Y(end); 
    d = (x-mup_)./ volp_ - theta*volp_; % Kappa in memoire
    grid(end,:,:) = (1-normcdf(d)).*c;

 %% on calcule les autre
 
for k=2:m
          
           % little trick to calculate the interval :
           x  = Y(k);
           d = (x-mup_)./ volp_ - theta*volp_;
           grid(k,:,:) = normcdf(d).*c;
           
           x = Y(k-1);
           d = (x-mup_)./ volp_ - theta*volp_; 
           grid(k,:,:) = squeeze(grid(k,:,:)) - normcdf(d) .* c;

end


%% check if grid is ok (check at 1,1)

% fun = @(x) exp(theta.*x) .* normpdf(x,mup_(1,1),volp_(1,1) );
% check = integral(fun,-Inf,Y(1)) ;
% diff = abs(check-grid(1,1,1));
% assert(diff < 1e-6  );







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     function  [grid] = create_grid3D_Mprime(m,reg,Y,theta,mup,phi,volp,M)

grid = zeros(m+1,m,reg);

%% add lag return dim to mu

       mup_ = (1-phi').*mup';
       lag = Y*phi';
       mup_ = lag+mup_;
       
       volp_ = repmat(volp',m,1);
       
      
%% compute constant c (in front of the intervals)

      c = exp(theta*mup_ + theta^2 * volp_.^2 /2) .* volp_;
      
%% function for derivative of normcdf N'

    N_prime = @(z) exp(-(z.^2)/2) / (sqrt(2*pi));
       
%% fill top and bottom borders (at q=0 and q=m)

    x  = Y(1);
    d = (x-mup_)./ volp_ - theta*volp_; % Kappa in memoire
    grid(1,:,:) = N_prime(d).*c;

    x  = Y(end); 
    d = (x-mup_)./ volp_ - theta*volp_; % Kappa in memoire
    grid(end,:,:) = N_prime(d).*c; % (1-N_prime(d)).*c;

 %% compute the others
 
for k=2:m
          
           % little trick to calculate the interval :
           x  = Y(k);
           d = (x-mup_)./ volp_ - theta*volp_;
           grid(k,:,:) = N_prime(d).*c;
           
           x = Y(k-1);
           d = (x-mup_)./ volp_ - theta*volp_; 
           grid(k,:,:) = squeeze(grid(k,:,:)) - N_prime(d) .* c;

end

%% add last piece

temp = repmat((mup_ + theta.*volp_.^2),1,1,m+1);
temp = permute(temp,[3 1 2]);
grid = temp.*M - grid;

%% check if grid is ok (check at 1,1 and 2,1)

% fun = @(x)  x .* exp(theta.*x) .* normpdf(x,mup_(1,1),volp_(1,1) );
% check1 = integral(fun,-Inf,Y(1)) ;
% check2 = integral(fun,Y(end),Inf) ;
% diff1 = abs(check1-grid(1,1,1));
% diff2 = abs(check2-grid(end,1,1));
% assert(diff1 < 1e-6  );
% assert(diff2 < 1e-6  );




    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     function  [grid] = create_grid5D_M(m,q,reg,S,Y,theta,mup,phi,volp)

grid = nan(m+1,q,m+1,q-1,reg); 

S1 = S(2:end);

%% add lag return dim to mu

       mup_ = (1-phi').*mup';
       lag = Y*phi';
       mup_ = lag+mup_;
       
       volp_ = repmat(volp',q,1);
       
%% cste in front of the integrals

       c = exp(theta*mup_ + theta^2 * volp_.^2 /2);
       
       % fill at s=0, s=0
       c_ = repmat(c,1,1,q-1);
       c_ = permute(c_,[1 3 2]);
       grid(1,:,1,:,:) = c_;  
       
%%


for i = 2:length(S)
           
        s_q = S(i); 
        x  = log(S1./s_q);
    
        for j = 1:length(Y)
            
            y_v = Y(j);
            mup = mup_(j,:);

            for k = 1:length(x)
                
                if k == 1 % lower border case
                    maxx = x(k);
                    Y_temp = min(Y,maxx);                
                else 
                    minx = x(k-1);
                    maxx = x(k);
                    Y_temp = max(Y,minx);
                    Y_temp = min(Y_temp,maxx);
                end
                                
                %%%%%% little trick to calculate the interval %%%%%
                   d = (Y_temp-mup)./ volp' - theta*volp';
                   temp = normcdf(d).*c(j,:);
                   temp = diff(temp);
                   grid(i,j,k,:,:) = temp ;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                               
            end
            
            % higher border case:
            k = m+1;
            minx = x(end);
            Y_temp = max(Y,minx);
            d = (Y_temp-mup)./ volp' - theta*volp';
            temp = (1-normcdf(d)).*c(j,:);
            temp = -diff(temp);
            grid(i,j,k,:,:) = temp ;
           
        end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     function  [grid] = create_grid5D_Mprime(m,q,reg,S,Y,theta,mup,phi,volp,M)

         
grid = nan(m+1,q,m+1,q-1,reg); 

S1 = S(2:end);

%% add lag return dim to mu

       mup_ = (1-phi').*mup';
       lag = Y*phi';
       mup_ = lag+mup_;
       
       volp_ = repmat(volp',q,1);
       
%% cste in front of the integrals

       c = exp(theta*mup_ + theta^2 * volp_.^2 /2) .* volp_;
       
       % fill at s=0, s=0
       c_ = repmat(c,1,1,q-1);
       c_ = permute(c_,[1 3 2]);
       grid(1,:,1,:,:) = c_; 

%% function for derivative of normcdf N'

    N_prime = @(z) exp(-(z.^2)/2) / (sqrt(2*pi));
    
%%

for i = 2:length(S)
           
        s_q = S(i); 
        x  = log(S1./s_q);
    
        for j = 1:length(Y)
            
            y_v = Y(j);
            mup = mup_(j,:);

            for k = 1:length(x)
                
                if k == 1 % lower border case
                    maxx = x(k);
                    Y_temp = min(Y,maxx);                
                else 
                    minx = x(k-1);
                    maxx = x(k);
                    Y_temp = max(Y,minx);
                    Y_temp = min(Y_temp,maxx);
                end
                                
                %%%%%%%%%%%%%%%
                % little trick to calculate the interval :
                   d = (Y_temp-mup)./ volp' - theta*volp';
                   temp = N_prime(d).*c(j,:);
                   temp = diff(temp);
                   grid(i,j,k,:,:) = temp ;
                %%%%%%%%%%%%%%%%
                               
            end
            
            % higher border case:
            k = m+1;
            if i ~= m+1
                minx = x(end);
                Y_temp = max(Y,minx);
                d = (Y_temp-mup)./ volp' - theta*volp';
                temp = N_prime(d).*c(j,:); %(1-N_prime(d)).*c(j,:);
                temp = -diff(temp);
                grid(i,j,k,:,:) = temp ;
            else
                minx = x(end);
                Y_temp = max(Y,minx);
                d = (Y_temp-mup)./ volp' - theta*volp';
                temp = N_prime(d).*c(j,:); %(1-N_prime(d)).*c(j,:);
                temp = diff(temp);
                grid(i,j,k,:,:) = temp ;
            end
           
        end
end

%% add last piece

temp = repmat((mup_ + theta.*volp_.^2),1,1,m+1,m+1,q-1);
temp = permute(temp,[3 1 4 5 2 ]); 
grid = temp.*M - grid;








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A,B] = create_A_B_2D(X,F,q,reg)
             
    A = zeros(q+1,reg);
    B = zeros(q+1,reg);
    
     dx = repmat((X(2:end)-X(1:end-1))',reg,1); 
     B(2:end-1,:) = ( F(2:end,:)-F(1:end-1,:) ) ./ dx';
     
     temp = repmat(X(1:(end-1))',reg,1);                 
     A(2:end-1,:) = F(1:end-1,:)- B(2:end-1,:).*temp'; 
       
     A(end,:) = A(end-1,:);
     B(end,:) = B(end-1,:);  
     A(1,:) = A(2,:);
     B(1,:) = B(2,:); 
     
     
     
     
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A,B1,B2,B3] = create_A_B123_3D(X,Y,F,m,q,reg)
             
    
      dx = repmat((X(2:end)-X(1:end-1))',reg,1); 
      dx_ = repmat(dx,1,1,q-1);
      dx_ = permute(dx_,[2 3 1]);
      B1 = ( F(2:end,1:end-1,:)-F(1:end-1,1:end-1,:) ) ./ dx_ ;
      
      dy = repmat((Y(2:end)-Y(1:end-1))',reg,1); 
      dy_ = repmat(dy,1,1,m);
      dy_ = permute(dy_,[3 2 1]);
      B2 = ( F(1:end-1,2:end,:)-F(1:end-1,1:end-1,:) ) ./ dy_ ;
      
      dxy = dx_ .* dy_ ;
      B3 = ( F(2:end,2:end,:) - F(2:end,1:end-1,:) - F(1:end-1,2:end,:) ...
                + F(1:end-1,1:end-1,:) ) ./dxy ;
      

      Y_ = repmat(Y(1:end-1,:),1,reg,m);
      Y_ = permute(Y_,[3 1 2]) ; 
      
      X_ = repmat(X(1:end-1),1,reg,q-1);
      X_ = permute(X_,[1 3 2]);     
      
      A = F(1:end-1,1:end-1,:) - B1.*X_ - B2.*Y_ + B3.*X_.*Y_ ;
 
      B1 = B1 - B3.*Y_;
     
      B2 = B2 - B3.*X_ ;     
      
      % add higher border
      A(m+1,:,:)  = A(m,:,:);
      B1(m+1,:,:) = B1(m,:,:);
      B2(m+1,:,:) = B2(m,:,:);
      B3(m+1,:,:) = B3(m,:,:);
      
          
     
     


   