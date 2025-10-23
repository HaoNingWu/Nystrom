
close all
clear all

% Spherical configurations and quadrature methods for integral equations of the second kind
% by C. An and H.-N. Wu

% MATLAB codes written by H.-N. Wu in Aug 2024
% Please add Chebfun and the sphere_approx_toolbox_v3.0 onto path

% Validation point set
Xt = get_Xt( ); Xt = Xt';



%% Equation to be solved
example_idx = 7;
% 1 3 5 7 for reproducing

% singulat kernel
% h_idx = 1; % non-singular kernel h = 1
% h_idx = 2: h(|x-y|) = |x-y|^{nu}
% h_idx = 3: h(|x-y|) = log(|x-y|)
% h_idx = 4: h(|x-y|) = |x-y|^{nu} |x+y|^{nu2}

switch example_idx
    case 1
        h_idx = 1; % non-singular kernels
        % specify nu as needed
        nu = 0; nu2 = 0; % place holding; not necessary to change
        Kernel = @(x,y) real(sin(10*norm(x-y,2)));
        f = @(i) sin(10*sqrt(2*(1-i))); I = integral(f,-1,1);
        func = @(x,y,z) 1- 2*pi*I;
        % Kernel = @(x,y) 1;
        % func = @(x,y,z) 1-4*pi; % testing for h = k =1;
    case 3
        h_idx = 2; % non-singular kernels
        % specify nu as needed
        nu = -.5; nu2 = 0;
        Kernel = @(x,y) real(cos(10*norm(x-y,2)));
        % f = @(i) sqrt(2*(1-i)).^(nu).*cos(10*sqrt(2*(1-i))); I = integral(f,-1,1);
        f = @(u) sqrt(u) .* cos(10*u);
        I = integral(f, 0, 2);
        % Kernel = @(x,y) 1;
        % f = @(i) sqrt(2*(1-i)).^(nu); I = integral(f,-1,1);
        func = @(x,y,z) 1- 2*pi*I;
        % func = @(x,y,z) 1 - pi*(4*log(2)-2); % testing for k = 1, h = log
    case 5
        h_idx = 3; % non-singular kernels
        % specify nu as needed
        nu = 0; nu2 = 0;
        Kernel = @(x,y) 1;
        % f = @(i) sqrt(2*(1-i)).^(nu).*sin(2*sqrt(2*(1-i))); I = integral(f,-1,1);
        func = @(x,y,z) 1 - pi*(4*log(2)-2); % testing for k = 1, h = log
        % n = 5;
        % t = floor(1.5*n);
        % m = (t+1)^2;
    case 7
        h_idx = 4; % non-singular kernels
        % specify nu as needed
        nu = -0.5; nu2 = -0.5;
        % Kernel = @(x,y) 1;
        % f = @(i) sqrt(2*(1-i)).^(nu).*sqrt(2*(1+i)).^(nu2); I = integral(f,-1,1);
        % f = @(i) sqrt(2*(1-i)).^(nu).*sqrt(2*(1+i)).^(nu2).*sin(10*sqrt(2*(1-i))); I = integral(f,-1,1);
        Kernel = @(x,y) real(sin(10*norm(x-y,2)));
        % Kernel = @(x,y) 1;
        % f = @(i) sqrt(2*(1-i)).^(nu).*sqrt(2*(1+i)).^(nu2).*sin(10*sqrt(2*(1-i))); I = integral(f,-1,1);

        f = @(theta) sin(20 * sin(theta/2)) .* sqrt(sin(theta/2) .* cos(theta/2));
        I = integral(f, 0, pi);

        func = @(x,y,z) 1- 2*pi*I;

end



n = 20;
k = 1;
Kmat = [];Kmatt = []; Error_max = [];
for t = 15:1:35
    % t = floor(1.2*n);
    m = (t+1)^2;
    ts(k) = t;
    for point_idx = 1:4

        switch point_idx
            case 1 % Spherical t-designs
                X_k = loadStd(t,(t+1)^2);
            case 2% Minimal energy points
                X_k = loadME( t, (t+1)^2 );
            case 3  % Fekete points
                X_k = loadMD( t, (t+1)^2 );
            case 4 % Equal area points
                X_k = sca_data_cap(2, 1 ,m, 0);
                X_k = X_k';


        end


        % RHS f sampling
        f = ones(m,1).*func(X_k(:,1),X_k(:,2),X_k(:,3));
        ft = ones(size(Xt,1),1)*func(Xt(:,1),Xt(:,2),Xt(:,3));



        %% Numerical Scheme Stage 1
        % Computing weights W_j
        Y_L = getQ( X_k, n )';
        W = Weights(n,nu,nu2,Y_L,Y_L,h_idx);

        % Kernel matrix K
        for i = 1:m
            for j = 1:m
                Kmat(i, j) = Kernel(X_k(i,:),X_k(j,:));
            end
        end

        % matrix in the linear system
        WKmat = W.*Kmat;

        % Solving the system
        phi = (eye(m) - WKmat) \ f;


        Cond1(k,point_idx) = cond(eye(m) - WKmat);


        %% Numerical Scheme Stage 2
        M = size(Xt,1);
        Yt_eqp = get_Yt( n, Xt' );
        Wt = Weights(n,nu,nu2,Y_L,Yt_eqp',h_idx);

        for i = 1:M
            for j = 1:m
                Kmatt(i, j) = Kernel(Xt(i,:),X_k(j,:));
            end
        end


        WKmatt = Wt.*Kmatt;
        Phi = WKmatt * phi + ft;
        Error_max(k,point_idx) = norm(Phi -  ones(size(Xt,1),1),inf);

    end
    fprintf('First subplot in progress: %.2f%%\n', (k/size(15:1:35,2)) * 100);
    k = k+1;

end

% Define markers, markersize, line width, and colors
markers = {'-o', '-s', '-^', '-d'};
markersize = 8;
linewidth = 2;
colors = {[227, 66, 52]/255, [0, 168, 107]/255, [0, 123, 167]/255, [255, 191, 0]/255};  % RGB colors

subplot(1,2,1)

% Plot with markers, markersize, lines, and colors
semilogy(ts, Error_max(:,1), markers{1}, 'MarkerSize', markersize, 'LineWidth', linewidth, 'MarkerFaceColor', colors{1}, 'Color', colors{1});
hold on
semilogy(ts, Error_max(:,2), markers{2}, 'MarkerSize', markersize, 'LineWidth', linewidth, 'MarkerFaceColor', colors{2}, 'Color', colors{2});
semilogy(ts, Error_max(:,3), markers{3}, 'MarkerSize', markersize, 'LineWidth', linewidth, 'MarkerFaceColor', colors{3}, 'Color', colors{3});
semilogy(ts, Error_max(:,4), markers{4}, 'MarkerSize', markersize, 'LineWidth', linewidth, 'MarkerFaceColor', colors{4}, 'Color', colors{4});

% Add legend
set(gca, 'fontsize', 15)
legend('spherical $t$-designs', 'minimal energy points', 'Fekete points', 'equal area points', 'Interpreter', 'latex', 'FontSize', 20);

title(['$n=', num2str(n), '$ with varying $m=(t+1)^2$'], 'Interpreter', 'latex', 'FontSize', 26);
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 24);
ylabel('uniform error', 'Interpreter', 'latex', 'FontSize', 24);
grid on, box on









k = 1;
Kmat = [];Kmatt = []; Error_max = [];
for n = 15:1:35
    t = floor(1.2*n);
    m = (t+1)^2;
    ns(k) = n;
    for point_idx = 1:4

        switch point_idx
            case 1 % Spherical t-designs
                X_k = loadStd(t,(t+1)^2);
            case 2% Minimal energy points
                X_k = loadME( t, (t+1)^2 );
            case 3  % Fekete points
                X_k = loadMD( t, (t+1)^2 );
            case 4 % Equal area points
                X_k = sca_data_cap(2, 1 ,m, 0);
                X_k = X_k';
        end


        % RHS f sampling
        f = ones(m,1).*func(X_k(:,1),X_k(:,2),X_k(:,3));
        ft = ones(size(Xt,1),1)*func(Xt(:,1),Xt(:,2),Xt(:,3));



        %% Numerical Scheme Stage 1
        % Computing weights W_j
        Y_L = getQ( X_k, n )';
        W = Weights(n,nu,nu2,Y_L,Y_L,h_idx);

        % Kernel matrix K
        for i = 1:m
            for j = 1:m
                Kmat(i, j) = Kernel(X_k(i,:),X_k(j,:));
            end
        end

        % matrix in the linear system
        WKmat = W.*Kmat;

        % Solving the system
        phi = (eye(m) - WKmat) \ f;

        Cond2(k,point_idx) = cond(eye(m) - WKmat);



        %% Numerical Scheme Stage 2
        M = size(Xt,1);
        Yt_eqp = get_Yt( n, Xt' );
        Wt = Weights(n,nu,nu2,Y_L,Yt_eqp',h_idx);

        for i = 1:M
            for j = 1:m
                Kmatt(i, j) = Kernel(Xt(i,:),X_k(j,:));
            end
        end


        WKmatt = Wt.*Kmatt;
        Phi = WKmatt * phi + ft;
        Error_max(k,point_idx) = norm(Phi -  ones(size(Xt,1),1),inf);

    end
    fprintf('Second subplot in progress: %.2f%%\n', (k/size(15:1:35,2))*100);
    k = k+1;

end

% Define markers, markersize, line width, and colors
markers = {'-o', '-s', '-^', '-d'};
markersize = 8;
linewidth = 2;
colors = {[227, 66, 52]/255, [0, 168, 107]/255, [0, 123, 167]/255, [255, 191, 0]/255};  % RGB colors

subplot(1,2,2)
% Plot with markers, markersize, lines, and colors
semilogy(ns, Error_max(:,1), markers{1}, 'MarkerSize', markersize, 'LineWidth', linewidth, 'MarkerFaceColor', colors{1}, 'Color', colors{1});
hold on
semilogy(ns, Error_max(:,2), markers{2}, 'MarkerSize', markersize, 'LineWidth', linewidth, 'MarkerFaceColor', colors{2}, 'Color', colors{2});
semilogy(ns, Error_max(:,3), markers{3}, 'MarkerSize', markersize, 'LineWidth', linewidth, 'MarkerFaceColor', colors{3}, 'Color', colors{3});
semilogy(ns, Error_max(:,4), markers{4}, 'MarkerSize', markersize, 'LineWidth', linewidth, 'MarkerFaceColor', colors{4}, 'Color', colors{4});

% Add legend
set(gca, 'fontsize', 15)
legend('spherical $t$-designs', 'minimal energy points', 'Fekete points', 'equal area points', 'Interpreter', 'latex', 'FontSize', 20);

title(['varying $n$ with $m = (\lfloor 1.2n\rfloor+1)^2$'], 'Interpreter', 'latex', 'FontSize', 26);
xlabel('degree $n$ of hyperinterpolants', 'Interpreter', 'latex', 'FontSize', 24);
ylabel('uniform error', 'Interpreter', 'latex', 'FontSize', 24);
grid on, box on






figure,

subplot(1,2,1)

% Plot with markers, markersize, lines, and colors
semilogy(ts, Cond1(:,1), markers{1}, 'MarkerSize', markersize, 'LineWidth', linewidth, 'MarkerFaceColor', colors{1}, 'Color', colors{1});
hold on
semilogy(ts, Cond1(:,2), markers{2}, 'MarkerSize', markersize, 'LineWidth', linewidth, 'MarkerFaceColor', colors{2}, 'Color', colors{2});
semilogy(ts, Cond1(:,3), markers{3}, 'MarkerSize', markersize, 'LineWidth', linewidth, 'MarkerFaceColor', colors{3}, 'Color', colors{3});
semilogy(ts, Cond1(:,4), markers{4}, 'MarkerSize', markersize, 'LineWidth', linewidth, 'MarkerFaceColor', colors{4}, 'Color', colors{4});
% Add legend
set(gca, 'fontsize', 15)
legend('spherical $t$-designs', 'minimal energy points', 'Fekete points', 'equal area points', 'Interpreter', 'latex', 'FontSize', 20);

title('$n=20$ with varying $m=(t+1)^2$', 'Interpreter', 'latex', 'FontSize', 26);
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 24);
ylabel('condition number', 'Interpreter', 'latex', 'FontSize', 24);
grid on, box on



subplot(1,2,2)
% Plot with markers, markersize, lines, and colors
plot(ns, Cond2(:,1), markers{1}, 'MarkerSize', markersize, 'LineWidth', linewidth, 'MarkerFaceColor', colors{1}, 'Color', colors{1});
hold on
plot(ns, Cond2(:,2), markers{2}, 'MarkerSize', markersize, 'LineWidth', linewidth, 'MarkerFaceColor', colors{2}, 'Color', colors{2});
plot(ns, Cond2(:,3), markers{3}, 'MarkerSize', markersize, 'LineWidth', linewidth, 'MarkerFaceColor', colors{3}, 'Color', colors{3});
plot(ns, Cond2(:,4), markers{4}, 'MarkerSize', markersize, 'LineWidth', linewidth, 'MarkerFaceColor', colors{4}, 'Color', colors{4});

% Add legend
set(gca, 'fontsize', 15)
legend('spherical $t$-designs', 'minimal energy points', 'Fekete points', 'equal area points', 'Interpreter', 'latex', 'FontSize', 20);

title(['varying $n$ with $m = (\lfloor 1.2n\rfloor+1)^2$'], 'Interpreter', 'latex', 'FontSize', 26);
xlabel('degree $n$ of hyperinterpolants', 'Interpreter', 'latex', 'FontSize', 24);
ylabel('condition number', 'Interpreter', 'latex', 'FontSize', 24);
grid on, box on

