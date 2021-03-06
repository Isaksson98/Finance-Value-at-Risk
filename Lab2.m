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

VaR_95 = Vp * norminv(0.95) * sigma_p;
VaR_975 = Vp * norminv(0.975) * sigma_p;
VaR_99 = Vp * norminv(0.99) * sigma_p;

%% 1 b)

r = zeros(size(u,1),1); %logaritmerade portföljavkastningar

for i = 1:size(u,1)
        r(i) = log(1 + weights' * u(i,:)');     
end

sigma_EWMA = EWMA(r, r(2)^2, 0.94);

VaR_t_95 = zeros(length(u),1);
VaR_t_99 = zeros(length(u),1);

VaR_t_95(502:length(u)) = (1 - exp(-norminv(0.95) * sigma_EWMA(502:length(u)))) * Vp;
VaR_t_99(502:length(u)) = (1 - exp(-norminv(0.99) * sigma_EWMA(502:length(u)))) * Vp;

figure(1)

plot(VaR_t_95)
title("konfidensnivå 95%")
figure(2)

plot(VaR_t_99)
title("konfidensnivå 99%")

%% 1 c)

percentile_95 = floor(500 * (1 - 0.95));
percentile_99 = floor(500 * (1 - 0.99));

VaR_95_rullande = zeros(length(u),1);
VaR_99_rullande = zeros(length(u),1);

delta_Vp = zeros(1, length(u));

for i = 501:length(u)
delta_Vp = u(i - 499:i,:) * (weights * Vp);
delta_VP_sorted = sort(delta_Vp);

VaR_95_rullande(i + 1) = delta_VP_sorted(percentile_95);
VaR_99_rullande(i + 1) = delta_VP_sorted(percentile_99);
end

delta_VP_sorted = sort(delta_Vp);
ES_95 = mean(delta_VP_sorted(1:25));

figure(3)
plot(-VaR_95_rullande)
title("konfidensnivå 95%")

figure(4)
plot(-VaR_99_rullande)
title("konfidensnivå 99%")

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
    
R_star = R .* sigma_i(1:length(sigma_i) - 1) / sigma_i(i + 1);
    
delta_Vp = R_star(i - 499:i,:) * Vp;
delta_VP_sorted = sort(delta_Vp);

VaR_95_hull(i + 1) = delta_VP_sorted(percentile_95); %ska man ta abs här?
VaR_99_hull(i + 1) = delta_VP_sorted(percentile_99);
end

figure(3)
plot(-VaR_95_hull)
title("konfidensnivå 95%")

figure(4)
plot(-VaR_99_hull)
title("konfidensnivå 99%")

%% 1 e)

%för b
LP_b = r * Vp;
[h0_b95_alpha025, ER_b95_alpha025] = calc_error_rate(LP_b(502:end), VaR_t_95(502:length(VaR_t_95)), 0.95, 0.05 / 2);
[h0_b95_alpha005, ER_b95_alpha005] = calc_error_rate(LP_b(502:end), VaR_t_95(502:length(VaR_t_95)), 0.95, 0.01 / 2);

[h0_b99_alpha025, ER_b99_alpha025] = calc_error_rate(LP_b(502:end), VaR_t_99(502:length(VaR_t_99)), 0.99, 0.05 / 2);
[h0_b99_alpha005, ER_b99_alpha005] = calc_error_rate(LP_b(502:end), VaR_t_99(502:length(VaR_t_99)), 0.99, 0.01 / 2);

%för c
LP_c = u * weights * Vp;
[h0_c95_alpha025, ER_c95_alpha025] = calc_error_rate(LP_c(502:end), -VaR_95_rullande(502:length(VaR_95_rullande)-1), 0.95, 0.05 / 2);
[h0_c95_alpha005, ER_c95_alpha005] = calc_error_rate(LP_c(502:end), -VaR_95_rullande(502:length(VaR_95_rullande)-1), 0.95, 0.01 / 2);

[h0_c99_alpha025, ER_c99_alpha025] = calc_error_rate(LP_c(502:end), -VaR_99_rullande(502:length(VaR_99_rullande)-1), 0.99, 0.05 / 2);
[h0_c99_alpha005, ER_c99_alpha005] = calc_error_rate(LP_c(502:end), -VaR_99_rullande(502:length(VaR_99_rullande)-1), 0.99, 0.01 / 2);

%för c
LP_d = u * weights * Vp;
[h0_d95_alpha025, ER_d95_alpha025] = calc_error_rate(LP_d(502:end), -VaR_95_hull(502:length(VaR_95_hull)-1), 0.95, 0.05 / 2);
[h0_d95_alpha005, ER_d95_alpha005] = calc_error_rate(LP_d(502:end), -VaR_95_hull(502:length(VaR_95_hull)-1), 0.95, 0.01 / 2);

[h0_d99_alpha025, ER_d99_alpha025] = calc_error_rate(LP_d(502:end), -VaR_99_hull(502:length(VaR_99_hull)-1), 0.99, 0.05 / 2);
[h0_d99_alpha005, ER_d99_alpha005] = calc_error_rate(LP_d(502:end), -VaR_99_hull(502:length(VaR_99_hull)-1), 0.99, 0.01 / 2);

%% 1 f)

%b)
test_storhet_b95 = serialdependance(VaR_t_95(502:length(VaR_t_95)),LP_b(502:end))
test_storhet_b99 = serialdependance(VaR_t_95(502:length(VaR_t_95)),LP_b(502:end))

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

problem_3 = xlsread(filename, 'Problem 3');

%Risk factors
SP500 = problem_3(:,1);
sigma_VIX = problem_3(:,2);
rate = problem_3(:,3);

valuation_date = '2-Feb-2021';
S = 3826.31;
r = 0.021;
div = 0.05;

K1 = 3800; %call
T1 = days252bus(valuation_date, '21-Mar-2021') / 252;
Bid1 = 20.80;
Ask1 = 20.99;
Holdings1 = 10000;

K2 = 3750; %put
T2 = days252bus(valuation_date, '21-Apr-2021') / 252;
Bid2 = 22.67;
Ask2 = 22.81;
Holdings2 = 10000;

K3 = 3850; %call
T3 = days252bus(valuation_date, '2021-Sep-21') / 252;
Bid3 = 21.88;
Ask3 = 22.03;
Holdings3 = 20000;

%Implicit_Sigma

implicit_sigma_1 = Implicit_Sigma(Ask1, S * exp(-div * T1), K1, r, T1, 1);
implicit_sigma_2 = Implicit_Sigma(Ask2, S * exp(-div * T2), K2, r, T2, 0);
implicit_sigma_3 = Implicit_Sigma(Ask3, S * exp(-div * T3), K3, r, T3, 1);

%Prices

Price1 = BSM(S * exp(-div * T1), K1, r, implicit_sigma_1, T1);
Price2 = BSM_put(S * exp(-div * T2), K2, r, implicit_sigma_2, T2);
Price3 = BSM(S * exp(-div * T3), K3, r, implicit_sigma_3, T3);

%Greeks

[delta1, vega1, rho1] = Greeks(S, K1, r, implicit_sigma_1, T1, div);
[delta2, vega2, rho2] = Greeks(S, K2, r, implicit_sigma_2, T2, div);
[delta3, vega3, rho3] = Greeks(S, K3, r, implicit_sigma_3, T3, div);

%
