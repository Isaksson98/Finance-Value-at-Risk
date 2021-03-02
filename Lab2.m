%% 1 a)
clear;
filename='timeSeries.xlsx';
portfolio = xlsread(filename, 'Problem 1 and 2');
portfolio = flip(portfolio);

u = zeros(size(portfolio,1), (size(portfolio,2)));

for i = 2:size(portfolio,1) %avkastningar
    for j = 1:size(portfolio,2)
        u(i,j) = (portfolio(i, j) - portfolio(i - 1, j)) / portfolio(i - 1, j);    
    end
end

cov_matrix = cov(u(2:size(u,1),:));

weights = ones(size(u,2), 1) / size(u,2);

sigma_p = sqrt(weights' * cov_matrix * weights);

Vp = 10^7;

VaR_95 = Vp * norminv(0.95) * sigma_p
VaR_975 = Vp * norminv(0.975) * sigma_p
VaR_99 = Vp * norminv(0.99) * sigma_p

%% 1 b)

r = zeros(size(u,1),1); %logaritmerade portföljavkastningar

for i = 1:size(u,1)
        r(i) = log(1 + weights' * u(i,:)');     
end

sigma_EWMA = EWMA(r, r(2)^2, 0.94);

VaR_t_95 = (1 - exp(-norminv(0.95) * sigma_EWMA(502:length(u)))) * Vp;
VaR_t_99 = (1 - exp(-norminv(0.99) * sigma_EWMA(502:length(u)))) * Vp;

%% 1 c)

percentile_95 = floor(500 * (1 - 0.95));
percentile_99 = floor(500 * (1 - 0.99));

VaR_95_rullande = zeros(length(u),1);
VaR_99_rullande = zeros(length(u),1);

delta_Vp = zeros(1, length(u));

for i = 501:length(u) %%fel värden - vet varför fel instoppat
delta_Vp = u(i - 499:i,:) * (weights * Vp);
delta_VP_sorted = sort(delta_Vp);

VaR_95_rullande(i + 1) = delta_VP_sorted(percentile_95); %konstga värden???
VaR_99_rullande(i + 1) = delta_VP_sorted(percentile_99);
end

delta_VP_sorted = sort(delta_Vp);
ES_95 = mean(delta_VP_sorted(1:25));

%% 1 d)

VaR_95_hull = zeros(length(u),1);
VaR_99_hull = zeros(length(u),1);
 
delta_Vp = zeros(1, 500);

R = u * weights;
R_mean = mean(R(2:21));

sum_i = 0;
for i = 2:21
    sum_i = sum_i + (R(i) - R_mean)^2;
end
sigma_init = sum_i / 19;

sigma_i = EWMA(R, sigma_init, 0.94);

    
for i = 501:length(R)
    
R_star = R .* sigma_i(1:length(sigma_i) - 1) / sigma_i(i);
    
delta_Vp = R_star(i - 499:i,:) * Vp;
delta_VP_sorted = sort(delta_Vp);

VaR_95_hull(i + 1) = delta_VP_sorted(percentile_95); %ska man ta abs här?
VaR_99_hull(i + 1) = delta_VP_sorted(percentile_99);
end

%% 1 e)

%för b
LP_b = r * Vp;
[h0_b95, error_rate_b95] = calc_error_rate(LP_b(502:end), VaR_t_95, 0.05);
[h0_b99, error_rate_b99] = calc_error_rate(LP_b(502:end), VaR_t_99, 0.01);

%för c
LP_c = u * weights * Vp;
[h0_c95, error_rate_c95] = calc_error_rate(LP_c(502:end), -VaR_95_rullande(502:length(VaR_95_rullande)-1), 0.05);
[h0_c99, error_rate_c99] = calc_error_rate(LP_c(502:end), -VaR_99_rullande(502:length(VaR_99_rullande)-1), 0.01);

%% 1 f)

%b)
test_storhet_b95 = serialdependance(VaR_t_95,LP_b(502:end))
test_storhet_b99 = serialdependance(VaR_t_99,LP_b(502:end))

%c)
test_storhet_c95 = serialdependance(-VaR_95_rullande(502:length(VaR_95_rullande)-1),LP_c(502:end))
test_storhet_c99 = serialdependance(-VaR_99_rullande(502:length(VaR_99_rullande)-1),LP_c(502:end))

%d)
test_storhet_d95 = serialdependance(-VaR_95_hull(502:length(VaR_95_hull)-1),LP_c(502:end))
test_storhet_d99 = serialdependance(-VaR_99_hull(502:length(VaR_99_hull)-1),LP_c(502:end))

%Critical limit
chi_95 = chi2inv(0.95,1);
chi_99 = chi2inv(0.99,1);


%% 2 a)

u_hat = prctile(sort(r), 95);
y = r-u_hat;
n_u = sum(r > u_hat)

[beta_xi, val] = fmincon(@(x)EVT(x,y, n_u), [0.9 0.4], [], [], [], [], [0 0], []);

beta = beta_xi(1);
xi = beta_xi(2);

VaR_EVT = -(u_hat + (beta / xi) * (((length(r) / n_u) * 0.95)^(-xi) - 1) ) * Vp

%b)

plot(sigma_EWMA)

average_volatility = zeros((length(sigma_EWMA) - 260), 1);

for i = 1:(length(average_volatility))
   average_volatility(i) = mean(sigma_EWMA(i:i + 260)); 
end
[argvalue, argmax] = max(average_volatility);
beginning = argmax;
ending = argmax + 260; %%hittar perioden (5år) med störst average volatilitet

u_hat = prctile(sort(r(beginning:ending)), 95);
y = r-u_hat;
n_u = sum(r > u_hat)

[beta_xi, val] = fmincon(@(x)EVT(x,y, n_u), [0.9 0.4], [], [], [], [], [0 0], []);

beta = beta_xi(1);
xi = beta_xi(2);

VaR_EVT_volatile = -(u_hat + (beta / xi) * (((length(r) / n_u) * 0.95)^(-xi) - 1) ) * Vp

%% 3 a)






