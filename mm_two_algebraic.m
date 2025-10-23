function mu = mm_two_algebraic(n,nu1,nu2)

x = chebfun('x');
for ell = 0:n
    p = legpoly(ell);
    mu_temp = 0;
    mu(ell+1) = 2^((nu1+nu2)/2)*2*pi * sum((1-x).^(nu1/2).*(1+x).^(nu2/2).*legpoly(ell),[-1,1]);
end

end
