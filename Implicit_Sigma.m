function sigma = Implicit_Sigma(ask_price, S, K, r, T, type)

max = 100;
min = 0;
middle = (max + min) / 2;
price = 0;

while(abs(price - ask_price) > 0.0001) %intervallhalvering'
    
    if type == 1
        price = BSM(S, K, r, middle, T);
    else
        price = BSM_put(S, K, r, middle, T);
    end
    
    if(price - ask_price > 0)
        max = middle;
    else
        min = middle;
    end
    
    middle = (max + min) / 2;
end

sigma = middle;
end

