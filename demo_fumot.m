% clc; clear;
N = 2^12 + 1; % disk, but requires sufficient regularity, larger the better.
theta = linspace(0, 2*pi, N);
theta = theta(1:(N-1));
nodes = [cos(theta); sin(theta)];
femm_opt = struct('deg', 3, 'qdeg', 8, 'min_area', 1e-3, 'edge', nodes);
gamma_opt = struct('X', 0.4, 'M', 0.6);
beta_opt  = struct('X', 0.15, 'M', 0.7, 'F', 0.2); 

% tau = gammaX / betaF, tau is not -1.
% mu  = betaX/betaF - 1.

opt = struct('femm_opt', femm_opt, 'reg', 1e-4, 'gamma', gamma_opt, 'beta', beta_opt);
fmt = FUMOT(opt);

p = struct('aF', fmt.parameter.aF);
[Q, u] = fmt.forward_ex(p,0);
fmt.backward_ex(Q, u)
