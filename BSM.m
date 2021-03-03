function C = BSM(S, K, r, sigma, T)
d1 = (log(S / K) + (r + (sigma^2) / 2) * T) / (sigma * T^0.5);
d2 = d1 - sigma * T^0.5;

C = S * normcdf(d1) - K * exp(-r * T) * normcdf(d2);
end