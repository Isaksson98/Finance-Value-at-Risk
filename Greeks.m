function [delta, vega, rho] = Greeks(S, K, r, sigma, T, q, call)

d1 = (log(S / K) + (r + (sigma^2) / 2) * T) / (sigma * T^0.5);
d2 = d1 - sigma * T^0.5;

if call == 1
    delta = exp(-q * T) * normcdf(d1);
    vega = S * exp(-q * T) * normpdf(d1) * sqrt(T);
    rho = K * T * exp(-r * T) * normcdf(d2);
else
    delta = -exp(-q * T) * normcdf(-d1);
    vega = S * exp(-q * T) * normpdf(d1) * sqrt(T);
    rho = -K * T * exp(-r * T) * normcdf(-d2);
end

end