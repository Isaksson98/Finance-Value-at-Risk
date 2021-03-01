function test_storhet = serialdependance(VaR, LP)

n00 = 0;
n01 = 0;
n10 = 0;
n11 = 0;

L = LP > VaR;

for i = 2:length(L)
    
    if L (i-1) == 0 && L(i) == 0
        n00 = n00+1;
    elseif  L (i-1) == 0 && L(i) == 1
        n01= n01+1;
    elseif L (i-1) == 1 && L(i) == 0
        n10 = n10+1;
    else 
        n11 = n11+1;
    end
end

pi01 = n01/(n01+n00);
pi11 = n11/(n10+n11);
pi_tot = (n01+n11)/(n01+n00+n10+n11);

test_storhet_1 = -2*log( ((1-pi_tot)^(n00+n10)) * pi_tot^(n01+n11));
test_storhet_2 = 2*log( ( (1-pi01)^n00) * (pi01^n01)*((1-pi11)^n10)*pi11^n11);
test_storhet = test_storhet_1+test_storhet_2;

end