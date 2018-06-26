function R = simulate_BS(mu,sigma,r,T,n)



R = (mu-.5*sigma^2-r)*(T/n) + randn(n,1) * sigma * sqrt(T/n);

