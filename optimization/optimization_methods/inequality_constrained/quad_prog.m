% An active-set nonbinding direction based method.
%
%   1. Solves the inequality constrained quadratic program
%
%           minimize g'x + (1/2)*x'*H*x
%           subject to Ax >= b
%
% 	2. If no initial subspace minimizer is provided, the program solves
%      a phase one linear program to determine an initial vertex.
%
%   3. Returns `unbounded' if the objective decreases in an unbounded
%      feasible direction.

% Jeremy Orosco
% MATH 271C, Spring 2016
% Coimbra Energy Group
% Department of Mechanical
% and Aerospace Engineering
% University of California, San Diego

% Copyright 2016 Jeremy Orosco
% This work is licensed under the
% Creative Commons Attribution-NonCommercial-ShareAlike 4.0
% International License. To view a copy of this license,
% visit http://creativecommons.org/licenses/by-nc-sa/4.0/.


function [x_star,y_star,r_star,f_star] = quad_prog(A,b,g,H,varargin)

    % start timer
    start_time = tic;
    
    % get system sizes
    n = size(A,2);
    m = size(A,1);
    
    % initialize / preallocate
    k = 0;                          % loop parameter
    tolerance = 10e01*sqrt(eps);    % zero radius
    iter_limit = 100;               % iteration limit
    alpha_h = zeros(m,1);           % feasible step storage
    not_updated = 1;                % p,q update status
    
    % get initial feasible point
    if nargin < 5
        x = phase_one(A,b);
        verbose = 1;
    else
        if numel(varargin) < 2
            x = phase_one(A,b,varargin{1});
            verbose = 1;
        else
            if numel(varargin{1}) == 0
                x = phase_one(A,b);
                verbose = varargin{2};
            else
                x = phase_one(A,b,varargin{1});
                verbose = varargin{2};
            end
        end
    end
    
    % compute residual
    r = A*x-b;
    
    % get active set
    active = find(abs(r) <= tolerance);
    Aa = A(active,:);
    mw = rank(Aa);
    
    % define initial working set
    [~,~,pa] = lu(Aa,'vector');
    working = active(pa(1:mw));
    Aw = A(working,:);
    
    % get minimum multiplier
    y = linsolve(Aw',g+H*x);
    [yks,s] = min(y);
    
    while (yks < -tolerance) && (k < iter_limit)
        
        % reinitialize constraint update indices
        removed = 0;
        added = 0;

        % output iteration data
        if verbose && k == 0
            fprintf('\n')
            fprintf('=============================================================\n')
            fprintf('  itn    min(y)      p''Hp         step    w-   w+      obj.  \n')
            fprintf('-------------------------------------------------------------\n')
            fprintf (' %3u   %9.2e       -           -       -    -   %9.2e\n', k, y(s),g'*x+(1/2)*x'*H*x);
        end
        
        % solve KKT system
        if not_updated
            es=zeros(n+mw,1);es(n+s)=1;
            K = [ H  Aw';
                 Aw  zeros(mw)];
            kr = es;
            pq = linsolve(K,kr);
            p = pq(1:n);
            q = -pq(n+1:n+mw);
        end
        
        % compute maximum feasible step
        if q(s) > tolerance
            alpha_s = -(y(s)/q(s));
        else
            alpha_s = inf;
        end
        for i = 1:m
           if A(i,:)*p < -tolerance
               alpha_h(i) = r(i)/(-A(i,:)*p);
           else
               alpha_h(i) = Inf;
           end
        end
        [alpha_f,t] = min(alpha_h);
        alpha = min(alpha_f,alpha_s);
        
        % check for unbounded solution
        if alpha == Inf
            disp('unbounded')
            break
        end
        
        % update parameters
        x = x + alpha*p;
        y = y + alpha*q;
        r = r + alpha*A*p;
        
        % check for independence of blocking constraint
        if alpha_f < alpha_s
            
            K = [ H Aw';
                 Aw zeros(mw)];
            kr = [A(t,:)';zeros(mw,1)];
            uv = linsolve(K,kr);
            u = uv(1:n);
            v = uv(n+1:end);
            
            % check blocking constraint dependence, update multipliers
            if norm(u) < tolerance || mw == n
                gamma = y(s)/v(s);
                y = [y-gamma*v;gamma];
                not_updated = 1;
            else
                rho = -(A(t,:)*p)/(A(t,:)*u);
                p = p + rho*u;
                q = [q - rho*v;rho];
                y = [y+0;0];
            end
          
            working = [working+0;t];
            mw = mw + 1;
            added = t;
            
        end
        if abs(y(s)) < tolerance
            
            ws = working(s);
            dws = find(working ~= ws);
            y = y(dws);
            working = working(dws);
            mw = mw - 1;
            removed = ws;
            
        end
        
        [~,s] = min(y);
        
        Aw = A(working,:);
        
        k = k + 1;
        
        if verbose && k <= 20
            fprintf (' %3u   %9.2e   %9.2e   %9.2e  %2u   %2u   %9.2e\n', k, y(s), p'*H*p, alpha, removed, added, g'*x+(1/2)*x'*H*x);
        end
        
        y_star = zeros(m,1);y_star(working)=y;yks = min(y_star);
    end
    
    % assign outputs
    x_star = x;
    y_star = zeros(m,1);y_star(working)=y;
    r_star = A*x_star-b;
    f_star = g'*x_star + (1/2)*x_star'*H*x_star;
    
    if verbose
        % display elapsed time
        fprintf('\n-> time elapsed: %0.4f\n',toc(start_time))

        % display iterations
        fprintf('\n-> method iterations: %u\n',k)

        % display optimizing value of x
        fprintf('\n-> optimizing value of x:\n\n')
        disp(round(x_star,5))

        % display optimum
        fprintf('-> objective optimum: %0.4f\n\n',f_star)
        fprintf('=============================================================\n\n')
    end

end



