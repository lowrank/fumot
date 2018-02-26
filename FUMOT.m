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
        
        
        function S = forward_em(obj, u0, noise)
            if nargin == 2
                noise = 0;
            end
            
            A = obj.cache.sm + obj.cache.mm;
            w = zeros(obj.cache.n, 1);
            l = obj.parameter.eta .* obj.parameter.aF .* u0; %load vector
            ql = obj.mapping(l, obj.model.space.elems, obj.model.facet.ref');
            tL = obj.model.build('l', ql);
            b =  tL;
            w(obj.cache.dof) = A(obj.cache.dof, obj.cache.dof) \ b(obj.cache.dof);
            
            [DX, DY] = obj.model.builder.reference_grad(obj.model.rffs.nodes);
            [wx, wy] = obj.model.gradient(w, DX, DY);
            % now ux, uy are close to the actual values. We randomly choose
            % one for the gradient.
            wgrad = zeros(obj.cache.n, 2);
            for i = 1:size(obj.model.space.elems, 2)
                wgrad(obj.model.space.elems(:, i), 1) = wx(:, i);
                wgrad(obj.model.space.elems(:, i), 2) = wy(:, i);
            end
            
            v = obj.load.em;
            v(obj.cache.dof) = 0;
            b = -A * v;
            v(obj.cache.dof) = A(obj.cache.dof, obj.cache.dof)\ b(obj.cache.dof);
            [vx, vy] = obj.model.gradient(v, DX, DY);
            vgrad = zeros(obj.cache.n, 2);
            for i = 1:size(obj.model.space.elems, 2)
                vgrad(obj.model.space.elems(:, i), 1) = vx(:, i);
                vgrad(obj.model.space.elems(:, i), 2) = vy(:, i);
            end
            
            C = obj.cache.sx + obj.cache.mx + obj.cache.mf;
            cl = obj.parameter.eta .* obj.parameter.aF .* v;
            qcl =  obj.mapping(cl, obj.model.space.elems, obj.model.facet.ref');
            b = obj.model.build('l', qcl);
            p = zeros(obj.cache.n, 1);
            p(obj.cache.dof) = C(obj.cache.dof, obj.cache.dof)\ b(obj.cache.dof);
            
            [px, py] = obj.model.gradient(p, DX, DY);
            pgrad = zeros(obj.cache.n, 2);
            for i = 1:size(obj.model.space.elems, 2)
                pgrad(obj.model.space.elems(:, i), 1) = px(:, i);
                pgrad(obj.model.space.elems(:, i), 2) = py(:, i);
            end
            
            [ux, uy] = obj.model.gradient(u0, DX, DY);
            % now ux, uy are close to the actual values. We randomly choose
            % one for the gradient.
            ugrad = zeros(obj.cache.n, 2);
            for i = 1:size(obj.model.space.elems, 2)
                ugrad(obj.model.space.elems(:, i), 1) = ux(:, i);
                ugrad(obj.model.space.elems(:, i), 2) = uy(:, i);
            end
            
            S = obj.parameter.gammaM * obj.parameter.dM .* (wgrad(:,1).*vgrad(:,1) + wgrad(:,2) .* vgrad(:, 2)) ...
                + obj.parameter.betaM * obj.parameter.aM .* w.*v ...
                - obj.parameter.betaF * obj.parameter.eta .* obj.parameter.aF .* u0 .* v ...
                +obj.parameter.gammaX * obj.parameter.dX .* (ugrad(:,1) .* pgrad(:,1) + ugrad(:,2) .* pgrad(:,2)) ...
                + (obj.parameter.betaX * obj.parameter.aX + obj.parameter.betaF * obj.parameter.aF) .* u0 .* p;
            
            S = S .*(1 + 0.5 * noise * (2 * rand(size(S)) - 1));
            
        end
        
        function S = forward_em_private(obj, maF, meta, u0)
            tic;
            A = obj.cache.sm + obj.cache.mm;
            w = zeros(obj.cache.n, 1);
            
            l = meta .* maF .* u0; %load vector
            
            ql = obj.mapping(l, obj.model.space.elems, obj.model.facet.ref');
            tL = obj.model.build('l', ql);
            b =  tL;
            w(obj.cache.dof) = A(obj.cache.dof, obj.cache.dof) \ b(obj.cache.dof);
            
            [DX, DY] = obj.model.builder.reference_grad(obj.model.rffs.nodes);
            [wx, wy] = obj.model.gradient(w, DX, DY);
            % now ux, uy are close to the actual values. We randomly choose
            % one for the gradient.
            wgrad = zeros(obj.cache.n, 2);
            for i = 1:size(obj.model.space.elems, 2)
                wgrad(obj.model.space.elems(:, i), 1) = wx(:, i);
                wgrad(obj.model.space.elems(:, i), 2) = wy(:, i);
            end
            
            
            v = obj.load.em;
            v(obj.cache.dof) = 0;
            b = -A * v;
            v(obj.cache.dof) = A(obj.cache.dof, obj.cache.dof)\ b(obj.cache.dof);
            [vx, vy] = obj.model.gradient(v, DX, DY);
            vgrad = zeros(obj.cache.n, 2);
            for i = 1:size(obj.model.space.elems, 2)
                vgrad(obj.model.space.elems(:, i), 1) = vx(:, i);
                vgrad(obj.model.space.elems(:, i), 2) = vy(:, i);
            end
            
            
            
            C = obj.cache.sx + obj.cache.mx + obj.cache.mf;
            cl = meta .* maF .* v;
            qcl =  obj.mapping(cl, obj.model.space.elems, obj.model.facet.ref');
            b = obj.model.build('l', qcl);
            p = zeros(obj.cache.n, 1);
            p(obj.cache.dof) = C(obj.cache.dof, obj.cache.dof)\ b(obj.cache.dof);
            
            [px, py] = obj.model.gradient(p, DX, DY);
            pgrad = zeros(obj.cache.n, 2);
            for i = 1:size(obj.model.space.elems, 2)
                pgrad(obj.model.space.elems(:, i), 1) = px(:, i);
                pgrad(obj.model.space.elems(:, i), 2) = py(:, i);
            end
            
            
            [ux, uy] = obj.model.gradient(u0, DX, DY);
            % now ux, uy are close to the actual values. We randomly choose
            % one for the gradient.
            ugrad = zeros(obj.cache.n, 2);
            for i = 1:size(obj.model.space.elems, 2)
                ugrad(obj.model.space.elems(:, i), 1) = ux(:, i);
                ugrad(obj.model.space.elems(:, i), 2) = uy(:, i);
            end
            
            S = obj.parameter.gammaM * obj.parameter.dM .* (wgrad(:,1).*vgrad(:,1) + wgrad(:,2) .* vgrad(:, 2)) ...
                + obj.parameter.betaM * obj.parameter.aM .* w.*v ...
                -obj.parameter.betaF * meta .* maF .* u0 .* v ...
                +obj.parameter.gammaX * obj.parameter.dX .* (ugrad(:,1) .* pgrad(:,1) + ugrad(:,2) .* pgrad(:,2)) ...
                + (obj.parameter.betaX * obj.parameter.aX + obj.parameter.betaF * maF) .* u0 .* p;
            toc;
        end
        
        function backward_ex_chk(obj, Q, u)
            assert(obj.parameter.betaF ~= 0); % special case I.
            tau = obj.parameter.gammaX / obj.parameter.betaF;
            mu = obj.parameter.betaX / obj.parameter.betaF - 1;
            assert(tau + 1 ~= 0); % special case II.

%             LapSqrtDx = LaplaceSqrtDiffusionFX(obj.model.space.nodes)';
%             lambda = LapSqrtDx ./ sqrt(obj.parameter.dX);
%             
%             NGradDx = NormedGradientDiffusionFX(obj.model.space.nodes)';
%             kappa  = (NGradDx ./ obj.parameter.dX).^2;
            
            
%             a = (obj.parameter.dX).^(-tau);
            a = (obj.parameter.dX);
%             b =  (1 + tau) * a .* (lambda - 0.25 * tau * kappa - ...
%                 mu * obj.parameter.aX./obj.parameter.dX);
            b = -(1 + tau) * mu * obj.parameter.aX;
            c =  (1 + tau) * Q./obj.parameter.betaF;
            
            % solve the nonlinear PDE. 
            % check first.
            tqdx = obj.mapping(a, obj.model.space.elems, obj.model.facet.ref');
            tqax = obj.mapping(b, obj.model.space.elems, obj.model.facet.ref');
            
            tmp = c.* u.^(-2);
            tmp2 = obj.mapping(tmp, obj.model.space.elems, obj.model.facet.ref');
            
            tS = obj.model.build('s', tqdx);
            tM = obj.model.build('m', tqax);
            tM2 = obj.model.build('m', tmp2);
            
            L = (u).^(tau + 1);
            L(obj.cache.dof) = 0;
            A = tS + tM + tM2;
            rhs = -A * L;
            L(obj.cache.dof) = A(obj.cache.dof, obj.cache.dof) \ rhs(obj.cache.dof);
            
            trisurf(obj.model.space.elems(1:3,:)', obj.model.space.nodes(1,:),...
                obj.model.space.nodes(2,:),L ./((u).^(tau + 1)));
            
            v = L ./(u).^(tau + 1);
            
            assert(norm(v-1) /norm(v) < 1e-7);      
        end
        
        function [aF, u] = backward_ex(obj, Q)
            % still needs theory to show this discretized nonlinear
            % equation is solvable and matches the solution. Depends on
            % what the paper is focusing on.
            assert(obj.parameter.betaF ~= 0); % special case I.
            tau = obj.parameter.gammaX / obj.parameter.betaF;
            mu = obj.parameter.betaX / obj.parameter.betaF - 1;
            assert(tau + 1 ~= 0); % special case II.

%             LapSqrtDx = LaplaceSqrtDiffusionFX(obj.model.space.nodes)';
%             lambda = LapSqrtDx ./ sqrt(obj.parameter.dX);
%             
%             NGradDx = NormedGradientDiffusionFX(obj.model.space.nodes)';
%             kappa  = (NGradDx ./ obj.parameter.dX).^2;
            
            
%             a = (obj.parameter.dX).^(-tau);
            a = (obj.parameter.dX);
%             b =  (1 + tau) * a .* (lambda - 0.25 * tau * kappa - ...
%                 mu * obj.parameter.aX./obj.parameter.dX);
            b = -(1 + tau) * mu * obj.parameter.aX;
            c =  (1 + tau) * Q./obj.parameter.betaF;
            
            tqdx = obj.mapping(a, obj.model.space.elems, obj.model.facet.ref');
            tqax = obj.mapping(b, obj.model.space.elems, obj.model.facet.ref');
            
            
            % test b and c.
            assert(all(b > 0));
            theta = (1 - tau) / (tau + 1);
            
            assert(all(theta * c <= 0));
            
            
            tS = obj.model.build('s', tqdx);
            tM = obj.model.build('m', tqax);
            
            tqc = obj.mapping(c, obj.model.space.elems, obj.model.facet.ref');
            tN = obj.model.build('m', tqc);
            
            % begin iteration. 
            % first, get initialized value for u.
            psi = (obj.load.ex).^(tau + 1); % boundary does not change.
            psi(obj.cache.dof) = 0;
            % the initial guess!!
            if tau >= 1
                A = tS + tM;
                rhs = -A * psi;
                psi(obj.cache.dof) = A(obj.cache.dof, obj.cache.dof) \ rhs(obj.cache.dof);  
            elseif tau < -1
                nu = max(psi);
                nu = -theta * nu^(-(1+theta));
                A = tS + tM + nu* tN;
                rhs = -A * psi;
                psi(obj.cache.dof) = A(obj.cache.dof, obj.cache.dof) \ rhs(obj.cache.dof);
            else
                A = tS + tM;
                rhs = -A * psi;
                psi(obj.cache.dof) = A(obj.cache.dof, obj.cache.dof) \ rhs(obj.cache.dof);  
                kappa = min(psi);
                nu = -theta * kappa^(-(1+theta));
            end

            
            err = 1e99;Iter = 0;
            while (err > 1e-12)  
                if (Iter == 0)
                    fprintf('\t Iteration \t error \n');
                end
                Iter = Iter + 1;

                if (tau >= 1)
                    tqcp = obj.mapping(c.*(psi).^(-(1+theta)), obj.model.space.elems, obj.model.facet.ref');
                    tL = obj.model.build('m', tqcp);
                    B = tS + tM + tL;
                    v = (obj.load.ex).^(tau + 1); % boundary does not change.
                    v(obj.cache.dof) = 0;
           
                    RHS = -B * v;
                    v(obj.cache.dof) = B(obj.cache.dof, obj.cache.dof)\ RHS(obj.cache.dof);
                elseif (tau < -1)
                    ld = obj.mapping(-c.*(psi.^(-theta) - nu *psi), obj.model.space.elems, obj.model.facet.ref');
                    LL = obj.model.build('l', ld);
                    B = tS + tM + nu * tN;
                    v = (obj.load.ex).^(tau + 1); % boundary does not change.
                    v(obj.cache.dof) = 0;
                    RHS = -B * v + LL;
                    v(obj.cache.dof) = B(obj.cache.dof, obj.cache.dof) \ RHS(obj.cache.dof);
                else
                    ld = obj.mapping(-c.*(psi.^(-theta) - nu *psi), obj.model.space.elems, obj.model.facet.ref');
                    LL = obj.model.build('l', ld);
                    
                    kappa = min(psi);
                    nu = -theta * kappa^(-(1+theta));
                    
                    B = tS + tM + nu * tN;
                    v = (obj.load.ex).^(tau + 1); % boundary does not change.
                    v(obj.cache.dof) = 0;
                    RHS = -B * v + LL;
                    v(obj.cache.dof) = B(obj.cache.dof, obj.cache.dof) \ RHS(obj.cache.dof);  
                end
                
                delta = psi - v;
                err = norm(delta)/norm(psi);
                psi = v;
                fprintf('\t %6d  \t %6.2e\n', Iter, err);

            end                
            
%             u = (psi).^(1/(tau+1))./(sqrt(obj.parameter.dX));
            u = (psi).^(1/(tau+1));
             
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
        
        function [res,flag,relres,iter,resvec] = backward_em(obj, S, maF, u0)
            ff = @(meta)(obj.forward_em_private(maF, meta, u0));
            [res,flag,relres,iter,resvec] = gmres(ff, S, 10, 1e-5, 5);
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

