function [h_0, error_rate] = calc_error_rate(L_P, VaR, p)

L = L_P > VaR;
X_T = sum(L);
error_rate = X_T / length(L);

T = length(L);

Z = (X_T - T * p) / sqrt(T * p * (1 - p));
h_0 = Z > norminv(1 - p);

end