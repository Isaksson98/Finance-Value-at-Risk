function f = EVT(x, y, n_u)

    beta = x(1);
    xi = x(2);

    f = 0;
    for i = 1:n_u
        f = f + log( (1/beta)*(1+xi*y(i)/beta)^(-1/xi-1) );
    end
    f = -f;
end