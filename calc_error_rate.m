function [h_0, error_rate] = calc_error_rate(L_P, VaR, c, alpha)

L = L_P < -VaR;
X_T = sum(L);
error_rate = X_T / length(L);

T = length(L);

p = 1 - c;

Z = (X_T - T * p) / sqrt(T * p * (1 - p));
h_0 = Z < norminv(1 - alpha / 2) && Z < norminv(alpha / 2);

end