function [delta, vega, rho] = Greeks(S, K, r, sigma, T, q)

d1 = (log(S / K) + (r + (sigma^2) / 2) * T) / (sigma * T^0.5);
d2 = d1 - sigma * T^0.5;

delta = exp(-q * T) * norm(d1);
vega = S * exp(-q * T) * norm(d1) * sqrt(T);
rho = K * T * exp(-r * T) * norm(d2);

end