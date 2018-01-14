% clc; clear;
N = 2^12 + 1; % disk, but requires sufficient regularity, larger the better.
theta = linspace(0, 2*pi, N);
theta = theta(1:(N-1)) + pi/(N-1);
nodes = [cos(theta); sin(theta)];
femm_opt = struct('deg', 2, 'qdeg', 6, 'min_area', 1e-3, 'edge', nodes);
gamma_opt = struct('X', 2, 'M', 0.6);
beta_opt  = struct('X', 0.15, 'M', 0.7, 'F', 0.25); 

% tau = gammaX / betaF, tau is not -1.
% mu  = betaX/betaF - 1.

opt = struct('femm_opt', femm_opt, 'reg', 1e-4, 'gamma', gamma_opt, 'beta', beta_opt);
fmt = FUMOT(opt);
%% add some noise.
tic;[Q, ex_sol] = fmt.forward_ex(0.0);toc;

%%
% fmt.backward_ex_chk(Q, ex_sol);

%% get absorption coefficient for fluoresence.
% caution, it is slow, 
tic;[aF, u] = fmt.backward_ex(Q);toc;
            
figure(3);
trisurf(fmt.model.space.elems(1:3,:)', fmt.model.space.nodes(1,:),...
                fmt.model.space.nodes(2,:),Q);

%%
fprintf('The relative error is %6.2e.\n', norm(aF - fmt.parameter.aF)/norm(fmt.parameter.aF));
figure(1);
trisurf(fmt.model.space.elems(1:3,:)', fmt.model.space.nodes(1,:),...
                fmt.model.space.nodes(2,:),aF);


figure(2);
trisurf(fmt.model.space.elems(1:3,:)', fmt.model.space.nodes(1,:),...
                fmt.model.space.nodes(2,:),u);


