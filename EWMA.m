function sigmas_p = EWMA(r, sigma_init, lambda)

n = size(r,1) + 1;

sigma_p_sq = zeros(n, 1);
sigma_p_sq(3) = sigma_init;

for i = 3:(n - 1)
    sigma_p_sq(i + 1) = (1 - lambda) * r(i)^2 + lambda * sigma_p_sq(i);
end

sigmas_p = sqrt(sigma_p_sq);
end