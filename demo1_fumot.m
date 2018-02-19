% clc; clear;
N = 2^2 + 1; % disk, but requires sufficient regularity, larger the better.
theta = linspace(0, 2*pi, N);
theta = theta(1:(N-1)) + pi/(N-1);
nodes = [cos(theta); sin(theta)];
femm_opt = struct('deg', 4, 'qdeg', 8, 'min_area', 1e-4, 'edge', nodes);
gamma_opt = struct('X', -2.6, 'M', -2.4);
beta_opt  = struct('X', -0.6, 'M', -0.4, 'F', -0.8); 

% tau = gammaX / betaF, tau is not -1.
% mu  = betaX/betaF - 1.

fprintf('tau is %6.2e.\n', gamma_opt.X / beta_opt.F);
fprintf('mu is %6.2e.\n', beta_opt.X / beta_opt.F - 1);

opt = struct('femm_opt', femm_opt, 'reg', 1e-4, 'gamma', gamma_opt, 'beta', beta_opt);
fmt = FUMOT(opt);
%% add some noise.
tic;[Q, u0] = fmt.forward_ex(0.02);toc;
tic;[S] = fmt.forward_em(u0, 0.02);toc;
% fmt.backward_ex_chk(Q, u0);

%% get absorption coefficient for fluoresence.
tic;[aF, u] = fmt.backward_ex(Q);toc;
figure(1);
fmt.plot(aF);
%% output of 1st stage
fprintf('The relative error is %6.2e.\n', norm(aF - fmt.parameter.aF , 1)/norm(fmt.parameter.aF, 1));
%% get quantum efficiency 
tic;[eta, flag,relres,iter,resvec] = fmt.backward_em(S, aF, u);toc;
%% output of 2nd stage
fprintf('The relative error is %6.2e.\n', norm(eta - fmt.parameter.eta)/norm(fmt.parameter.eta));
figure(2);
fmt.plot(eta);


