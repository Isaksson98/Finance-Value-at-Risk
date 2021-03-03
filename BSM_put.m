function P = BSM_put(S, K, r, sigma, T)
d1 = (log(S / K) + (r + (sigma^2) / 2) * T) / (sigma * T^0.5);
d2 = d1 - sigma * T^0.5;

P = K * exp(-r * T) * normcdf(-d2) - S * normcdf(-d1);
end