classdef FUMOT < handle
% Fluoresecence Ultrasound Modulated Optical Tomography
    
    properties (Access = public)
        cache
        measurement
        source
        load
        parameter
        model
        reg
    end
    
    methods
        function obj = FUMOT(opt)
            assert(isfield(opt, 'femm_opt'));
            assert(isfield(opt, 'gamma'));
            assert(isfield(opt, 'beta'));
            
            obj.cache = struct('s',[], 'm', [], 'n', [], 'dof', [], 'ndof',[],...
                'sx', [], 'mx', [], 'sm', [], 'mm', [], 'mf',[]); % cached variables
            obj.parameter = struct('dX', [], 'aX', [], 'dM', [], 'aM', [], ...
                'aF', [], 'eta', [], 'gammaX', 0, 'gammaM', 0, 'betaX', 0, ...
                'betaM', 0, 'betaF', 0);
            obj.source = struct('ex',[], 'em',[]); % emission source for auxillary function use
            obj.load = struct('ex', [], 'em',[]); 
            
            obj.model = femm(opt.femm_opt);
            
            obj.cache.s = obj.model.build('s', 1);
            obj.cache.m = obj.model.build('m', 1);
            
            obj.cache.n = size(obj.model.space.nodes, 2);
            obj.cache.ndof = unique(obj.model.space.edges);
            obj.cache.dof = setdiff(1:obj.cache.n, obj.cache.ndof);
            
            obj.reg = opt.reg;
            
            obj.parameter.dX = diffusionFX(obj.model.space.nodes)';
            obj.parameter.dM = diffusionFM(obj.model.space.nodes)';
            obj.parameter.aX = absorptionFX(obj.model.space.nodes)';
            obj.parameter.aM = absorptionFM(obj.model.space.nodes)';
            obj.parameter.aF = absorptionFF(obj.model.space.nodes)'; % to be recovered.
            obj.parameter.eta = quantumF(obj.model.space.nodes)'; % to be recovered.
            
            obj.parameter.gammaX = opt.gamma.X;
            obj.parameter.gammaM = opt.gamma.M;
            obj.parameter.betaX  = opt.beta.X;
            obj.parameter.betaM  = opt.beta.M;
            obj.parameter.betaF  = opt.beta.F;
            
            qdX = obj.mapping(obj.parameter.dX, obj.model.space.elems, obj.model.facet.ref');
            qdM = obj.mapping(obj.parameter.dM, obj.model.space.elems, obj.model.facet.ref');
            qaX = obj.mapping(obj.parameter.aX, obj.model.space.elems, obj.model.facet.ref');
            qaM = obj.mapping(obj.parameter.aM, obj.model.space.elems, obj.model.facet.ref');
            qaF = obj.mapping(obj.parameter.aF, obj.model.space.elems, obj.model.facet.ref');
            
            obj.cache.sx = obj.model.build('s', qdX);
            obj.cache.sm = obj.model.build('s', qdM);
            obj.cache.mx = obj.model.build('m', qaX);
            obj.cache.mm = obj.model.build('m', qaM);
            obj.cache.mf = obj.model.build('m', qaF);
           

            
            obj.source.ex = ex_source(obj.model.space.nodes); % full information
            obj.source.em = em_source(obj.model.space.nodes); % full information
            
            obj.load.ex = zeros(obj.cache.n, 1); % init
            obj.load.em = zeros(obj.cache.n, 1); % init
            
            obj.load.ex(obj.cache.ndof) = obj.source.ex(obj.cache.ndof); % boundary
            obj.load.em(obj.cache.ndof) = obj.source.em(obj.cache.ndof); % boundary
            
        end
        
        function [Q , u] = forward_ex(obj, noise)
            if nargin == 1
                noise = 0;
            end
            
            A = obj.cache.sx + obj.cache.mx + obj.cache.mf;
            u = obj.load.ex;
            b = -A * u;
            u(obj.cache.dof) = A(obj.cache.dof, obj.cache.dof) \ b(obj.cache.dof);
            
            [DX, DY] = obj.model.builder.reference_grad(obj.model.rffs.nodes);
            [ux, uy] = obj.model.gradient(u, DX, DY);
            % now ux, uy are close to the actual values. We randomly choose
            % one for the gradient.
            grad = zeros(obj.cache.n, 2);
            for i = 1:size(obj.model.space.elems, 2)
                grad(obj.model.space.elems(:, i), 1) = ux(:, i);
                grad(obj.model.space.elems(:, i), 2) = uy(:, i);
            end
            Q = obj.parameter.gammaX * obj.parameter.dX .* (grad(:,2).^2 + grad(:,1).^2) +...
                (obj.parameter.betaX * obj.parameter.aX +...
                obj.parameter.betaF * obj.parameter.aF) .* (u.^2);
            
            Q = Q .* (1 + 0.5 * noise * (2 * rand(size(Q)) - 1)); % add noise

        end
        
        function backward_ex_chk(obj, Q, u)
            assert(obj.parameter.betaF ~= 0); % special case I.
            tau = obj.parameter.gammaX / obj.parameter.betaF;
            mu = obj.parameter.betaX / obj.parameter.betaF - 1;
            assert(tau + 1 ~= 0); % special case II.

            LapSqrtDx = LaplaceSqrtDiffusionFX(obj.model.space.nodes)';
            lambda = LapSqrtDx ./ sqrt(obj.parameter.dX);
            
            NGradDx = NormedGradientDiffusionFX(obj.model.space.nodes)';
            kappa  = (NGradDx ./ obj.parameter.dX).^2;
            
            
            a = (obj.parameter.dX).^(-tau);
            b =  (1 + tau) * a .* (lambda - 0.25 * tau * kappa - ...
                mu * obj.parameter.aX./obj.parameter.dX);
            c =  (1 + tau) * a .* Q./obj.parameter.betaF;
            
            % solve the nonlinear PDE. 
            % check first.
            tqdx = obj.mapping(a, obj.model.space.elems, obj.model.facet.ref');
            tqax = obj.mapping(b, obj.model.space.elems, obj.model.facet.ref');
            
            tmp = c.* u.^(-2) ./ obj.parameter.dX;
            tmp2 = obj.mapping(tmp, obj.model.space.elems, obj.model.facet.ref');
            
            tS = obj.model.build('s', tqdx);
            tM = obj.model.build('m', tqax);
            tM2 = obj.model.build('m', tmp2);
            
            L = (sqrt(obj.parameter.dX) .* u).^(tau + 1);
            L(obj.cache.dof) = 0;
            A = tS + tM + tM2;
            rhs = -A * L;
            L(obj.cache.dof) = A(obj.cache.dof, obj.cache.dof) \ rhs(obj.cache.dof);
            
            trisurf(obj.model.space.elems(1:3,:)', obj.model.space.nodes(1,:),...
                obj.model.space.nodes(2,:),L ./((sqrt(obj.parameter.dX) .* u).^(tau + 1)));
            
            v = L ./((sqrt(obj.parameter.dX) .* u).^(tau + 1));
            
            assert(norm(v-1) /norm(v) < 1e-7);        
        end
        
        function [aF, psi] = backward_ex(obj, Q)
            % still needs theory to show this discretized nonlinear
            % equation is solvable and matches the solution. Depends on
            % what the paper is focusing on.
            assert(obj.parameter.betaF ~= 0); % special case I.
            tau = obj.parameter.gammaX / obj.parameter.betaF;
            mu = obj.parameter.betaX / obj.parameter.betaF - 1;
            assert(tau + 1 ~= 0); % special case II.

            LapSqrtDx = LaplaceSqrtDiffusionFX(obj.model.space.nodes)';
            lambda = LapSqrtDx ./ sqrt(obj.parameter.dX);
            
            NGradDx = NormedGradientDiffusionFX(obj.model.space.nodes)';
            kappa  = (NGradDx ./ obj.parameter.dX).^2;
            
            
            a = (obj.parameter.dX).^(-tau);
            b =  (1 + tau) * a .* (lambda - 0.25 * tau * kappa - ...
                mu * obj.parameter.aX./obj.parameter.dX);
            c =  (1 + tau) * a .* Q./obj.parameter.betaF;
            
            tqdx = obj.mapping(a, obj.model.space.elems, obj.model.facet.ref');
            tqax = obj.mapping(b, obj.model.space.elems, obj.model.facet.ref');
            
            
            % test b and c.
            assert(all(b > 0));
            theta = (1 - tau) / (tau + 1);
            
            assert(all(theta * c <= 0));
            
            
            tS = obj.model.build('s', tqdx);
            tM = obj.model.build('m', tqax);
            
            % begin iteration. 
            % first, get initialized value for u.
            psi = (sqrt(obj.parameter.dX) .* obj.load.ex).^(tau + 1);
            psi(obj.cache.dof) = 0;
            A = tS + tM;
            rhs = -A * psi;
            psi(obj.cache.dof) = A(obj.cache.dof, obj.cache.dof) \ rhs(obj.cache.dof);
            
            alpha = 0.99;
            nu = 0.1;
            kappa = 0.9;
            err = 1e99;Iter = 0;
            while (err > 1e-12)  
                if (Iter == 0)
                    fprintf('\t Iteration \t stepSize \t error \n');
                end
                Iter = Iter + 1;
                tmp2 = obj.mapping(c, obj.model.space.elems, obj.model.facet.ref');
                tM2 = obj.model.build('m', tmp2);
                
                A = tS + tM + sparse(1:obj.cache.n, 1:obj.cache.n, -theta * (psi).^(-theta-1)) * tM2;
                F = @(X)(( tS + tM ) * X + tM2 * (X).^(-theta));
                f = F(psi);
                delta = A(obj.cache.dof, obj.cache.dof) \ f(obj.cache.dof);
                     
                Ind = (delta > 0);
                z = psi(obj.cache.dof);
                if sum(Ind) 
                    lambda = min(alpha * z(Ind)./delta(Ind));
                    lambda = min(lambda, 1);
                else
                    lambda = 1;
                end
                
                psi(obj.cache.dof) = psi(obj.cache.dof) - lambda * delta;
                
                err = norm(delta)/norm(psi(obj.cache.dof));
                fprintf('\t %6d \t %6.2e\t %6.2e\n', Iter,lambda, err);

            end                
            
            u = (psi).^(1/(tau+1))./(sqrt(obj.parameter.dX));
             
            [DX, DY] = obj.model.builder.reference_grad(obj.model.rffs.nodes);
            [ux, uy] = obj.model.gradient(u, DX, DY);
            % now ux, uy are close to the actual values. We randomly choose
            % one for the gradient.
            grad = zeros(obj.cache.n, 2);
            for i = 1:size(obj.model.space.elems, 2)
                grad(obj.model.space.elems(:, i), 1) = ux(:, i);
                grad(obj.model.space.elems(:, i), 2) = uy(:, i);
            end
            aF = ((Q - obj.parameter.gammaX * obj.parameter.dX .* (grad(:,2).^2 + grad(:,1).^2))./(u.^2) -...
                (obj.parameter.betaX * obj.parameter.aX))./obj.parameter.betaF;
          
            % regularization ??
        end
        
        function plot(obj, f)
            % interpolation.
            N = 500;
            [X,Y] =meshgrid(linspace(-0.5,0.5, N));
            X = X(:);
            Y = Y(:);
            F = scatteredInterpolant(obj.model.space.nodes(1,:)', ...
                obj.model.space.nodes(2,:)', f);
            G = F(X, Y);
            surface(reshape(X, N,N), reshape(Y,N,N), reshape(G,N,N));
            shading interp; view(2);colorbar;colormap jet;
        end
    end
    
    methods(Static)
        function [interpolate] = mapping(func, elems, trans_ref)
            numberofqnodes = size(trans_ref, 1);
            interpolate = zeros(numberofqnodes, size(elems, 2));
            for i = 1: size(elems, 2)
                interpolate(:, i) = trans_ref * func(elems(:, i));
            end
        end
        function r = normsq(v)
            r = sum(v.^2);
        end

end
    
end

