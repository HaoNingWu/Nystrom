function mu = mm_log(n)
x = chebfun('x');
mu(1) = pi*(4*log(2)-2);
for ell = 1:n
    % p = legpoly(ell);
    mu_temp = 0;
    for k = 1:ell
        mu_temp = mu_temp + sum(x.^k.*legpoly(ell))/k;
    end
    mu(ell+1) = -pi * mu_temp;
    % mu(ell+1) = pi * sum(log(1-x).*legpoly(ell),[-1,1]);
end
end
