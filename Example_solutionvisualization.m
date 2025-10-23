
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
        % Kernel = @(x,y) real(sin(10*sqrt(2*(1-dot(x,y)))));
        Kernel = @(x,y) real(sin(10*norm(x-y,2)));
        f = @(i) sin(10*sqrt(2*(1-i))); I = integral(f,-1,1);
        func = @(x,y,z) 1- 2*pi*I;
        n = 20;
        t = floor(2*n);
        m = (t+1)^2;
    case 2
        h_idx = 1; % non-singular kernels
        % specify nu as needed
        nu = 0; nu2 = 0; % place holding; not necessary to change
        % Kernel = @(x,y) real(sin(10*sqrt(2*(1-dot(x,y)))));
        Kernel = @(x,y) real(sin(10*norm(x-y,2)));
        func = @(x,y,z) 0.75*exp(-((9*x-2).^2)/4-((9*y-2).^2)/4-((9*z-2).^2)/4) ...
            +0.75*exp(-((9*x+1).^2)/49-((9*y+1))/10-((9*z+1))/10)...
            +0.5*exp(-((9*x-7).^2)/4-((9*y-3).^2)/4-((9*z-5).^2)/4)...
            -0.2*exp(-((9*x-4).^2)-((9*y-7).^2)-((9*z-5).^2));
        n = 20;
        t = floor(2*n);
        m = (t+1)^2;
    case 3
        h_idx = 2; % non-singular kernels
        % specify nu as needed
        nu = -.5; nu2 = 0;
        Kernel = @(x,y) real(cos(10*norm(x-y,2)));
        % f = @(i) sqrt(2*(1-i)).^(nu).*cos(10*sqrt(2*(1-i))); I = integral(f,-1,1);

        f = @(u) sqrt(u) .* cos(10*u);
        I = integral(f, 0, 2);
        func = @(x,y,z) 1- 2*pi*I;
        n = 20;
        t = floor(2*n);
        m = (t+1)^2;
    case 4
        h_idx = 2; % non-singular kernels
        % specify nu as needed
        nu = -.9; nu2 = 0;
        Kernel = @(x,y) real(sin(10*norm(x-y,2)));
        % Kernel = @(x,y) 1;
        % f = @(i) sqrt(2*(1-i)).^(nu).*sin(10*sqrt(2*(1-i))); I = integral(f,-1,1);
        func = @(x,y,z) 0.75*exp(-((9*x-2).^2)/4-((9*y-2).^2)/4-((9*z-2).^2)/4) ...
            +0.75*exp(-((9*x+1).^2)/49-((9*y+1))/10-((9*z+1))/10)...
            +0.5*exp(-((9*x-7).^2)/4-((9*y-3).^2)/4-((9*z-5).^2)/4)...
            -0.2*exp(-((9*x-4).^2)-((9*y-7).^2)-((9*z-5).^2));
        n = 20;
        t = floor(2*n);
        m = (t+1)^2;
    case 5
        h_idx = 3; % non-singular kernels
        % specify nu as needed
        nu = 0; nu2 = 0;
        Kernel = @(x,y) 1;
        % f = @(i) sqrt(2*(1-i)).^(nu).*sin(2*sqrt(2*(1-i))); I = integral(f,-1,1);
        func = @(x,y,z) 1 - pi*(4*log(2)-2); % testing for k = 1, h = log
        n = 5;
        t = floor(2*n);
        m = (t+1)^2;
    case 6
        h_idx = 3; % non-singular kernels
        % specify nu as needed
        nu = 0; nu2 = 0;
        Kernel = @(x,y) 1;
        % f = @(i) sqrt(2*(1-i)).^(nu).*sin(2*sqrt(2*(1-i))); I = integral(f,-1,1);
        func = @(x,y,z) 0.75*exp(-((9*x-2).^2)/4-((9*y-2).^2)/4-((9*z-2).^2)/4) ...
            +0.75*exp(-((9*x+1).^2)/49-((9*y+1))/10-((9*z+1))/10)...
            +0.5*exp(-((9*x-7).^2)/4-((9*y-3).^2)/4-((9*z-5).^2)/4)...
            -0.2*exp(-((9*x-4).^2)-((9*y-7).^2)-((9*z-5).^2));
        n = 10;
        t = floor(2*n);
        m = (t+1)^2;
    case 7
        h_idx = 4; % non-singular kernels
        % specify nu as needed
        nu = -.5; nu2 = -0.5;
        Kernel = @(x,y) real(sin(10*norm(x-y,2)));
        % Kernel = @(x,y) 1;
        % f = @(i) sqrt(2*(1-i)).^(nu).*sqrt(2*(1+i)).^(nu2).*sin(10*sqrt(2*(1-i))); I = integral(f,-1,1);
        f = @(theta) sin(20 * sin(theta/2)) .* sqrt(sin(theta/2) .* cos(theta/2));
        I = integral(f, 0, pi);
        % f = @(i) sqrt(2*(1-i)).^(nu).*sqrt(2*(1+i)).^n (nu2); I = integral(f,-1,1);
        func = @(x,y,z) 1- 2*pi*I;
        n = 20;
        t = floor(2*n);
        m = (t+1)^2;
    case 8
        h_idx = 4; % non-singular kernels
        % specify nu as needed
        nu = -.5; nu2 = -0.5;
        % Kernel = @(x,y) real(sin(2*norm(x-y,2)));
        Kernel = @(x,y) 1;
        % f = @(i) sqrt(2*(1-i)).^(nu).*sqrt(2*(1+i)).^(nu2).*sin(2*sqrt(2*(1-i))); I = integral(f,-1,1);
        % f = @(i) sqrt(2*(1-i)).^(nu).*sqrt(2*(1+i)).^(nu2); I = integral(f,-1,1);
        func = @(x,y,z) 0.75*exp(-((9*x-2).^2)/4-((9*y-2).^2)/4-((9*z-2).^2)/4) ...
            +0.75*exp(-((9*x+1).^2)/49-((9*y+1))/10-((9*z+1))/10)...
            +0.5*exp(-((9*x-7).^2)/4-((9*y-3).^2)/4-((9*z-5).^2)/4)...
            -0.2*exp(-((9*x-4).^2)-((9*y-7).^2)-((9*z-5).^2));
        n = 10;
        t = floor(2*n);
        m = (t+1)^2;
end


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
    switch mod(example_idx,2)
        case 1
            f = ones(m,1).*func(X_k(:,1),X_k(:,2),X_k(:,3));
            ft = ones(size(Xt,1),1)*func(Xt(:,1),Xt(:,2),Xt(:,3));
        case 0
            ft = func(Xt(:,1),Xt(:,2),Xt(:,3));
            f=func(X_k(:,1),X_k(:,2),X_k(:,3));
    end



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

    % scale the plotting for better visual comparison
    x = Xt(:,1); y = Xt(:,2); z = Xt(:,3);
    tri = convhull([x y z]);

    C1 = Phi;
    C2 = abs(Phi - ones(size(Xt,1),1));


    switch mod(example_idx,2)
        case 1
            subplot(2,4,point_idx)
            [Fmax, imax] = max(C1);
            [Fmin, imin] = min(C1);
            scale = 0.1;
            FS = 1 + (scale/(Fmax-Fmin+eps))*(C1-Fmin);
            fg = trisurf(tri,x.*FS, y.*FS, z.*FS, C1,'facecolor','interp');
            set(fg,'EdgeColor', 'none');
            axis vis3d, axis equal tight,
            axis off
            set(gca, 'fontsize', 16)
            colorbar('eastoutside')
            switch point_idx
                case 1
                    title('spherical $t$-designs','interpreter','latex','fontsize', 24)
                case 2
                    title('minimal energy points','interpreter','latex','fontsize', 24)
                case 3
                    title('Fekete points','interpreter','latex','fontsize', 24)
                case 4
                    title('equal area points','interpreter','latex','fontsize', 24)
            end
            % colormap parula
            subplot(2,4,point_idx+4)
            [Fmax, imax] = max(C2);
            [Fmin, imin] = min(C2);
            scale = 0.7;
            FS = 1 + (scale/(Fmax-Fmin+eps))*(C2-Fmin);
            fg = trisurf(tri,x.*FS, y.*FS, z.*FS, C2,'facecolor','interp');
            set(fg,'EdgeColor', 'none');
            axis vis3d, axis equal tight,
            axis off
            set(gca, 'fontsize', 16)
            colorbar('eastoutside')
            title('absolute error','interpreter','latex','fontsize', 24)
            % colormap parula

        case 0
            subplot(1,4,point_idx)
            [Fmax, imax] = max(C1);
            [Fmin, imin] = min(C1);
            scale = 0.1;
            FS = 1 + (scale/(Fmax-Fmin+eps))*(C1-Fmin);
            fg = trisurf(tri,x.*FS, y.*FS, z.*FS, C1,'facecolor','interp');
            set(fg,'EdgeColor', 'none');
            axis vis3d, axis equal tight,
            axis off
            set(gca, 'fontsize', 16)
            colorbar('eastoutside')
            switch point_idx
                case 1
                    title('spherical $t$-designs','interpreter','latex','fontsize', 24)
                case 2
                    title('minimal energy points','interpreter','latex','fontsize', 24)
                case 3
                    title('Fekete points','interpreter','latex','fontsize', 24)
                case 4
                    title('equal area points','interpreter','latex','fontsize', 24)
            end
    end
end

