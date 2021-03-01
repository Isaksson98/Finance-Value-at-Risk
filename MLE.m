function num = MLE(x, y)

    beta = x(1);
    xi = x(2);

    n = floor(0.05*length(y));
    num = 0;
    for i = 1:n
        num = num + log( (1/beta)*(1+xi*y(i)/beta)^(-1/xi-1) );
    end
end