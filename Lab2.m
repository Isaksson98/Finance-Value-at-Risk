%% 1 a)
clear;
filename='timeSeries.xlsx';
portfolio = xlsread(filename, 'Problem 1 and 2', 'C:Q');
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

VaR_t_95(502:length(u)) = (1 - exp(-norminv(0.95) * sigma_EWMA(502:length(u))));
VaR_t_99(502:length(u)) = (1 - exp(-norminv(0.99) * sigma_EWMA(502:length(u))));

figure(1)
plot(VaR_t_95)
title("konfidensnivå 95%")
xlabel("t")
ylabel("VaR")

figure(2)
plot(VaR_t_99)
title("konfidensnivå 99%")
xlabel("t")
ylabel("VaR")
%% 1 c)

percentile_95 = floor(500 * (1 - 0.95));
percentile_99 = floor(500 * (1 - 0.99));

VaR_95_rullande = zeros(length(u),1);
VaR_99_rullande = zeros(length(u),1);

delta_Vp = zeros(1, length(u));

for i = 501:length(u)
delta_Vp = u(i - 499:i,:) * weights;
delta_VP_sorted = sort(delta_Vp);

VaR_95_rullande(i + 1) = -delta_VP_sorted(percentile_95);
VaR_99_rullande(i + 1) = -delta_VP_sorted(percentile_99);
end

delta_VP_sorted = sort(delta_Vp);
ES_95 = mean(delta_VP_sorted(1:25));

figure(3)
plot(VaR_95_rullande)
title("konfidensnivå 95%")
xlabel("t")
ylabel("VaR")

figure(4)
plot(VaR_99_rullande)
title("konfidensnivå 99%")
xlabel("t")
ylabel("VaR")

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
    
delta_Vp = R_star(i - 499:i,:);
delta_VP_sorted = sort(delta_Vp);

VaR_95_hull(i + 1) = -delta_VP_sorted(percentile_95);
VaR_99_hull(i + 1) = -delta_VP_sorted(percentile_99);
end

figure(5)
plot(VaR_95_hull)
title("konfidensnivå 95%")
xlabel("t")
ylabel("VaR")

figure(6)
plot(VaR_99_hull)
title("konfidensnivå 99%")
xlabel("t")
ylabel("VaR")
%% 1 e)

%för b
LP_b = u * weights;
[h0_b95_alpha025, ER_b95_alpha025] = calc_error_rate(LP_b(502:end), VaR_t_95(502:length(VaR_t_95)), 0.95, 0.05 / 2);
[h0_b95_alpha005, ER_b95_alpha005] = calc_error_rate(LP_b(502:end), VaR_t_95(502:length(VaR_t_95)), 0.95, 0.01 / 2);

[h0_b99_alpha025, ER_b99_alpha025] = calc_error_rate(LP_b(502:end), VaR_t_99(502:length(VaR_t_99)), 0.99, 0.05 / 2);
[h0_b99_alpha005, ER_b99_alpha005] = calc_error_rate(LP_b(502:end), VaR_t_99(502:length(VaR_t_99)), 0.99, 0.01 / 2);

%för c
LP_c = u * weights;
[h0_c95_alpha025, ER_c95_alpha025] = calc_error_rate(LP_c(502:end), VaR_95_rullande(502:length(VaR_95_rullande)-1), 0.95, 0.05 / 2);
[h0_c95_alpha005, ER_c95_alpha005] = calc_error_rate(LP_c(502:end), VaR_95_rullande(502:length(VaR_95_rullande)-1), 0.95, 0.01 / 2);

[h0_c99_alpha025, ER_c99_alpha025] = calc_error_rate(LP_c(502:end), VaR_99_rullande(502:length(VaR_99_rullande)-1), 0.99, 0.05 / 2);
[h0_c99_alpha005, ER_c99_alpha005] = calc_error_rate(LP_c(502:end), VaR_99_rullande(502:length(VaR_99_rullande)-1), 0.99, 0.01 / 2);

%för d
LP_d = u * weights;
[h0_d95_alpha025, ER_d95_alpha025] = calc_error_rate(LP_d(502:end), VaR_95_hull(502:length(VaR_95_hull)-1), 0.95, 0.05 / 2);
[h0_d95_alpha005, ER_d95_alpha005] = calc_error_rate(LP_d(502:end), VaR_95_hull(502:length(VaR_95_hull)-1), 0.95, 0.01 / 2);

[h0_d99_alpha025, ER_d99_alpha025] = calc_error_rate(LP_d(502:end), VaR_99_hull(502:length(VaR_99_hull)-1), 0.99, 0.05 / 2);
[h0_d99_alpha005, ER_d99_alpha005] = calc_error_rate(LP_d(502:end), VaR_99_hull(502:length(VaR_99_hull)-1), 0.99, 0.01 / 2);

%% 1 f)

%b)
test_storhet_b95 = serialdependance(VaR_t_95(502:length(VaR_t_95)),LP_b(502:end));
test_storhet_b99 = serialdependance(VaR_t_99(502:length(VaR_t_99)),LP_b(502:end));

%c)
test_storhet_c95 = serialdependance(VaR_95_rullande(502:length(VaR_95_rullande)-1),LP_c(502:end));
test_storhet_c99 = serialdependance(VaR_99_rullande(502:length(VaR_99_rullande)-1),LP_c(502:end));

%d)
test_storhet_d95 = serialdependance(-VaR_95_hull(502:length(VaR_95_hull)-1),LP_c(502:end));
test_storhet_d99 = serialdependance(-VaR_99_hull(502:length(VaR_99_hull)-1),LP_c(502:end));

%Critical limit
chi_95 = chi2inv(0.95,1);
chi_99 = chi2inv(0.99,1);

%% 2 a)

sorted_losses = sort(-u * weights);
u_hat = prctile(sorted_losses, 95);

find_larger = sorted_losses(sorted_losses > u_hat);
y = find_larger - u_hat;
n_u = length(y);

[beta_xi, val] = fmincon(@(x)EVT(x,y, n_u), [0.9 0.4], [], [], [], [], [0 0], []);

beta = beta_xi(1);
xi = beta_xi(2);

VaR_EVT = u_hat + (beta / xi) * ((((length(u) - 1) / n_u) * (1 - 0.99))^(-xi) - 1);

%b)

average_volatility = zeros((length(sigma_EWMA) - 260), 1);

for i = 1:(length(average_volatility))
   average_volatility(i) = mean(sigma_EWMA(i:i + 260)); 
end
[argvalue, argmax] = max(average_volatility);
beginning = argmax;
ending = argmax + 260; %%hittar perioden (5år) med störst average volatilitet

sorted_losses_2 = sort(-u(beginning:ending,:) * weights);
u_hat_2 = prctile(sorted_losses_2, 95);

find_larger_2 = sorted_losses_2(sorted_losses_2 > u_hat_2);
y_2 = find_larger_2 - u_hat_2;
n_u_2 = length(y_2);

[beta_xi_2, val_2] = fmincon(@(x)EVT(x,y_2, n_u_2), [0.9 0.4], [], [], [], [], [0 0], []);

beta_2 = beta_xi_2(1);
xi_2 = beta_xi_2(2);

VaR_EVT_volatile = u_hat_2 + (beta_2 / xi_2) * ((((length(u(beginning:ending)) - 1) / n_u_2) * (1 - 0.99))^(-xi_2) - 1);

figure(7)
subplot(2,1,1)
plot(y, gppdf(y, xi, beta));
title 'beta = 0.0199, xi = 0.1476';
subplot(2,1,2)
plot(y_2, gppdf(y_2, xi_2, beta_2));
title 'beta = 0.0194, xi = 0.3373';
%% 3 a)

problem_3 = xlsread(filename, 'Problem 3', 'C:E');

%Risk factors
SP500 = flip(problem_3(:,1));
VIX = flip(problem_3(:,2));
libor = flip(problem_3(:,3));

valuation_date = '2-Feb-2021';
S = 3826.31;
rf = 0.021;
div = 0.05;

K1 = 3800; %call
T1 = days252bus(valuation_date, '21-Mar-2021') / 252;
implicit_volatility_bid1 = 20.80 / 100;
implicit_volatility_ask1 = 20.99 / 100;
Holdings1 = 10000;

K2 = 3750; %put
T2 = days252bus(valuation_date, '21-Apr-2021') / 252;
implicit_volatility_bid2 = 22.67 / 100;
implicit_volatility_ask2 = 22.81 / 100;
Holdings2 = 10000;

K3 = 3850; %call
T3 = days252bus(valuation_date, '21-Sep-2021') / 252;
implicit_volatility_bid3 = 21.88 / 100;
implicit_volatility_ask3 = 22.03 / 100;
Holdings3 = 20000;

%Prices

price1 = BSM(S * exp(-div * T1), K1, rf, implicit_volatility_ask1, T1);
price2 = BSM_put(S * exp(-div * T2), K2, rf, implicit_volatility_ask2, T2);
price3 = BSM(S * exp(-div * T3), K3, rf, implicit_volatility_ask3, T3);

%Greeks

[delta1, vega1, rho1] = Greeks(S, K1, rf, implicit_volatility_ask1, T1, div, 1);
[delta2, vega2, rho2] = Greeks(S, K2, rf, implicit_volatility_ask2, T2, div, 0);
[delta3, vega3, rho3] = Greeks(S, K3, rf, implicit_volatility_ask3, T3, div, 1);

%Stora G

G = [[delta1 * S, delta2 * S, delta3 * S]; [vega1, vega2, vega3]; [rho1, rho2, rho3]];

%H-streck

H_streck = 0;

%lambda

lambda = zeros(length(SP500), 3);

for i = 2:length(SP500)
    lambda(i,1) = log(SP500(i) / SP500(i - 1));
    lambda(i,2) = (VIX(i) - VIX(i - 1)) / 100;
    lambda(i,3) = (libor(i) - libor(i - 1)) / 100;
end

%kovariansmatris lambda

cov_lambda = cov(lambda(2:length(lambda),:));

%portföljvolatilitet

h = [10000, 10000, 20000];
V = price1 * h(1) + price2 * h(2) + price3 * h(3);

sigma_p_optioner = sqrt((1 / (V ^ 2)) * h * G' * cov_lambda * G * h');

% VaR

VaR_optioner_1d_99 = V * norminv(0.99) * sigma_p_optioner;

%% 3 b)

gh = (G' * cov_lambda * G * h') / (sigma_p_optioner * V);
gh_VaR = norminv(0.99) * (G' * cov_lambda * G * h') / (sigma_p_optioner * V);
test = gh_VaR' * h';

hf = G * h';

sigma_p_hf = sqrt((1 / V ^ 2) * hf' * cov_lambda * hf);

ghf_VaR = norminv(0.99) * (cov_lambda * hf) / (sigma_p_hf * V);

test2 = ghf_VaR' * hf;