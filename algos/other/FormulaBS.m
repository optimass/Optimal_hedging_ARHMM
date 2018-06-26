function [Call, Put] = FormulaBS( S, K, r, tau, sigma, yield)
% Evaluation of the call and put values using Black-Scholes formula
%
% Input 
%        S: Current price of the underlying asset;
%        K: Strike price; 
%        r: (continuous) short-time interest rate for the period;
%        tau: time to maturity;
%        sigma: volatility in the B-S model;
%        yield: continuous dividend rate.
%
% N.B. All parameters can be vectors. If there is more than one, then they
% must have the same dimension.
%
% Output
%         Call: value of the European call,
%         Put : value of the European put.
%
% N.B.: The time and rate must be on the same scale as the one used for
%       estimating sigma. For example, if sigma is estimated with daily
%       prices, then tau should be in days, and r should be the interest
%       rate for 1 day.

if nargin < 6
    yield = 0;
end

d1 = (log(S./K)+ (r-yield).*tau + (sigma.^2).* tau / 2)./ ( sigma.* sqrt( tau ) );
d2 = d1 - sigma.* sqrt( tau );

Call    =  exp( -yield.* tau ).*S.* normcdf( d1 ) -  exp( -r.* tau ).* K.* normcdf( d2 );
Put     =  Call - S.* exp( -yield.* tau ) + K.* exp( -r.* tau ); % Put-Call parity