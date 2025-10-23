function mu = mm_algebraic(n,nu)

% 
for ell = 0:n
    mu(ell+1) = 2^(nu+2)*pi*pochhammer(-nu/2,ell)*gamma((nu+2)/2)/gamma(ell+nu/2+2);
end

% x = chebfun('x');
% for ell = 0:n
%     p = legpoly(ell);
%     mu(ell+1) = 2^(nu/2)*2*pi * sum((1-x).^(nu/2).*legpoly(ell),[-1,1]);
% end

end
