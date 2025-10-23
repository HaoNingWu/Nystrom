function W = Weights(n,nu,nu2,Y_L, Y_val, h_idx)

m = size(Y_L,2);

% Y_L: column - basis on the same points

% W:
% raw: quadrature points varying
% col：given quad points, validation points varying


switch h_idx
    case 1
        W =  4*pi*ones(size(Y_val,2),size(Y_L,2))/m;
    case 2 % h(x) = |ξ − x|^{\nu}
        mu = mm_algebraic(n,nu);
        mms = Moments(n,Y_val,mu);
        W = 4*pi/m*mms'*Y_L;
    case 3 % h(x) = log|ξ − x|
        mu = mm_log(n);
        mms = Moments(n,Y_val,mu);
        W = 4*pi/m*mms'*Y_L;
    case 4 % h(x) = |ξ − x|^{\nu1}|ξ + x|^{\nu2}
        mu = mm_two_algebraic(n,nu,nu2);
        mms = Moments(n,Y_val,mu);
        W = 4*pi/m*mms'*Y_L;
end

end
