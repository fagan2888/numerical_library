% Newton-Lagrange method using an Armijo line-search
%
%   1. Determines update directions for objective arguments, x,
%      and Lagrange multipliers, y, by solving the KKT system.
%
%   2. Performs Armijo line-search for acceptable objective reduction.
%
%   3. Calls objective function, which must return
%
%       [
%        f(unction),
%        c(onstraint vector),
%        g(radient of objective),
%        J(acobian matrix of constraint gradients),
%        H(essian of Lagrangian, w.r.t. x),
%        G(radient of Lagrangian, w.r.t. (x,y))
%       ]
%
%   4. Uses gradient of Lagrangian function to determine convergence.

% Jeremy Orosco
% MATH 271B, Winter 2016
% Coimbra Energy Group
% Department of Mechanical
% and Aerospace Engineering
% University of California, San Diego

% Copyright 2016 Jeremy Orosco
% This work is licensed under the
% Creative Commons Attribution-NonCommercial-ShareAlike 4.0
% International License. To view a copy of this license,
% visit http://creativecommons.org/licenses/by-nc-sa/4.0/.


function [x_new,f_new,g_new] = newton_lagrange(f_name,x,y)

    % define optimization parameters
    nf = 0; % function evaluation counter
    k = 0; % method iteration counter
    gamma_c = 1/2;
    eta_s = 1/4;
    
    % define limits
    tolerance = 10e02*sqrt(eps);
    iter_limit = 100;
    
    % initialize func. values, sizes, and multipliers
    if nargin < 3
        [f,c,g,J,H,G,~] = feval(f_name,x);nf=nf+1;
        n = numel(x);
        m = numel(c);
        y = zeros(m,1);
    else
        [f,c,g,J,H,G,~] = feval(f_name,x,y);nf=nf+1;
        n = numel(x);
        m = numel(y);
    end
    
    % start timer
    start_time = tic;
    
    % run modified Newton-Lagrange method
    while norm(G) >= tolerance && k < iter_limit
        
        % solve the KKT system
        Kl = [H J'; J zeros(m)];
        Kr = -[g-J'*y; c];
        mk = linsolve(Kl,Kr);
        
        % define update directions
        pk = mk(1:n);
        qk = -mk(n+1:n+m);
        
        % initialize line search increment
        j = 0;

        % initialize alpha
        alpha_k = 1;
        
        % get initial candidate argument/multiplier update
        x_new = x + alpha_k*pk;
        y_new = y + alpha_k*qk;
        
        % get updated values
        [f_new,c_new,g_new,J_new,H_new,G_new,~] = feval(f_name,x_new,y_new);nf=nf+1;

        % perform Armijo line search
        while (norm(G_new) > (1-alpha_k*eta_s)*norm(G)) && (j < iter_limit/2)

            alpha_k = gamma_c*alpha_k;
            x_new = x + alpha_k*pk;
            y_new = y + alpha_k*qk;
            [f_new,c_new,g_new,J_new,H_new,G_new,~] = feval(f_name,x_new,y_new);nf=nf+1;
            j = j + 1;

        end

        % update values
        x = x_new;
        y = y_new;
        f = f_new;
        c = c_new;
        g = g_new;
        J = J_new;
        H = H_new;
        G = G_new;

        % update counters
        k = k+1;

    end
    
    % display elapsed time
    disp(toc(start_time))
    
    % display optimization parameters
    space = sprintf('\n');
    fprintf('\nnumber of method iterations:\n')
    disp(space)
    disp(k)
    fprintf('\nnumber of function evaluations:\n')
    disp(space)
    disp(nf)
    fprintf('\noptimizing value of x:\n')
    disp(space)
    disp(x)
    fprintf('\noptimizing value of y:\n')
    disp(space)
    disp(y)
    fprintf('\nobjective optimum:\n')
    disp(space)
    disp(f)
    fprintf('\nLagrangian gradient at optimum:\n')
    disp(space)
    disp(G)

end

