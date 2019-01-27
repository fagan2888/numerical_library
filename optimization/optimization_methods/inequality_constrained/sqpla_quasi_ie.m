% Multivariate (in primal/dual/slack variables) augmented Lagrangian
% SQP merit-function method with quasi-Newton approximation to Hessian
% (using BFGS update) and backtracking line-search.
%
%   1. Solves inequality constrained optimization problem
%
%               minimize f(x)
%             subject to c(x) >= 0
%
%      with a multivariate augmented Lagrangian serving as a merit-function
%      for an SQP method.
%
%   2. Requires only 1st-derivative (function and constraint) information.
%
%   3. Performs backtracking line-search for acceptable merit-function
%      reduction.
%
%   4. Calls objective function, which must return
%
%       [
%        f(unction),
%        c(onstraint vector),
%        g(radient of objective),
%        J(acobian matrix of constraint gradients)
%       ]
%
%   5. The variable 'sol' can be used to select method of KKT solution:
%
%       'w'   : general method -> J[Z Y] = [0 W], Y'Z=0
%       'tq'  : triangular factorization -> J[Z Y] = [0 T], Y'Z=0
%       'pbs' : basis factorization -> J[Z Y] = [B S], Y'Z=0
%
%   6. BFGS update includes conditions for sufficient positive curvature.
%      If possible, curvature is made sufficiently positve. Otherwise,
%      update is foregone.
%
%   7. Uses gradient of Lagrangian function to determine convergence.

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


function [x,f,g] = sqpla_quasi_ie(f_name,x,y,sol)

    % define optimization parameters
    nf = 0;                 % init. func. eval. counter
    k = 0;                  % init. method eval. counter
    l = 0;
    gamma_c = 1/2;          % line-search contraction factor
    eta_s = 1/4;            % sufficient decrease ratio
    rho_p = 0;              % merit/model-func. penalty parameter
    eta_h = 1/2;

    % start timer
    start_time = tic;

    % get initial values
    [f,c,g,J,~,~,~] = feval(f_name,x);nf=nf+1;

    % get system sizes
    n = numel(x);
    m = numel(c);

    % initialize multiplier estimates and slack variables
    if (nargin < 4) || numel(y) == 0
        if (nargin < 3) || numel(y) == 0
            y = 10e-03*ones(m,1);
            sol = 'tq';
        elseif nargin < 4
            sol = 'tq';
        end
    end
    s = 10e-03*ones(m,1);

    % define limits
    tolerance = 10e02*sqrt(eps);
    iter_limit = 100;

    % construct gradient and Hessian
    Gl = g - J'*y;
    B = eye(n);

    % run La-SQP method w/ quasi-Newton BFGS routine
    while (norm(Gl) >= tolerance) && (k < iter_limit)
        
        % output iteration parameters
        cv = c(c < -tolerance);
        if k == 0
            fprintf('\n')
            fprintf(' itn   nf    step    violation    optimal     nz      rho_p\n')
            fprintf('-----------------------------------------------------------\n')
            str1 = sprintf (' %3g %4g       - ', k, nf);
            if norm(c) < tolerance && norm(Gl) < tolerance
                str2 = sprintf ('   (%7.1e)   (%7.1e)    %2u', norm(cv), norm(Gl), l);
            elseif norm(c) < tolerance && norm(Gl) >= tolerance
                str2 = sprintf ('   (%7.1e)    %7.1e     %2u', norm(cv), norm(Gl), l);
            elseif norm(c) >= tolerance && norm(Gl) < tolerance
                str2 = sprintf ('    %7.1e   (%7.1e)     %2u', norm(cv), norm(Gl), l);
            elseif norm(c) >= tolerance && norm(Gl) >= tolerance
                str2 = sprintf ('    %7.1e     %7.1e     %2u', norm(cv), norm(Gl), l);
            end
            disp([str1 str2])
        end
        
        % solve QP subproblem
        pk = quad_prog(J,-c,g,B,[],0);

        % get working sets
        ch = c + J*pk;
        working = ch < tolerance;
        Jw = J(working,:);
        cw = c(working);
        
        % get descent direction normal
        Y = orth(Jw');
        pn = Y*linsolve(Jw*Y,-cw);
        
        % solve KKT system
        [~,yt,~] = solveKKT_quasi(g,B,cw,Jw,sol);
        y_new = zeros(m,1);y_new(working) = yt;
        
        % get parameter updates
%         y_new = linsolve(J',g+B*pk);
        qk = y_new - y;
        rk = J*pk+c-s;
        
        % update penalty parameter
        dphi = (1/2)*(pk'*B*pk)+2*qk'*(c-s)-rk'*y_new;
        if dphi > 0;
            rho_h = 0;
        else
            rho_h = (rk'*y_new-2*qk'*(c-s))/((c-s)'*(c-s));
        end
        rho_p = max(rho_h+tolerance,4*rho_p);

        % construct merit function
        fp = f-y'*(c-s)+(rho_p/2)*((c-s)'*(c-s));

        % initialize alpha
        alpha_k = 1;

        % get initial candidate argument update
        x_new = x+alpha_k*pk;
        s_new = s+alpha_k*rk;
        [f_new,c_new,g_new,J_new,~,~,~] = feval(f_name,x_new);nf=nf+1;

        % construct updated merit function
        fp_new = f_new-y_new'*(c_new-s_new)+(rho_p/2)*((c_new-s_new)'*(c_new-s_new));

        % initialize line-search increment
        j = 0;

        % perform line-search
        while fp-fp_new < eta_s*alpha_k*(pk'*B*pk-rk'*y_new+(rho_p*(c-s)+2*qk)'*(c-s)) && j < iter_limit/2
            alpha_k = gamma_c*alpha_k;
            x_new = x+alpha_k*pk;
            y_new = y+alpha_k*qk;
            s_new = s+alpha_k*rk;
            [f_new,c_new,g_new,J_new,~,~,~] = feval(f_name,x_new);nf=nf+1;
            fp_new = f_new-y_new'*(c_new-s_new)+(rho_p/2)*((c_new-s_new)'*(c_new-s_new));
            j = j + 1;
        end
        Gl_new = g_new-J_new'*y_new;

        % compute FDH update
        dk = x_new-x;
        wk = Gl_new-Gl;

        % compute BFGS update
        sigma = alpha_k*(1-eta_h)*pk'*B*pk+tolerance;
        if wk'*dk > sigma
            B = B+(wk*wk')/(wk'*dk)-(B*(dk*dk')*B')/(dk'*B*dk);
        else
            z = x+alpha_k*pn;
            dk = x_new-z;
            [~,~,g_z,J_z,~,~,~] = feval(f_name,z);nf=nf+1;
            Gl_z = g_z-J_z'*y_new;
            wk = Gl_new-Gl_z;
            if wk'*dk >= sigma
                B = B+(wk*wk')/(wk'*dk)-(B*(dk*dk')*B')/(dk'*B*dk);
            end
            l = l+1;
        end
        
        % update values
        x = x_new;
        y = y_new;
        f = f_new;
        c = c_new;
        g = g_new;
        J = J_new;
        
        % update slack variables
        if rho_p == 0;
            s = max(0,c);
        else
            s = max(0,c-y/rho_p);
        end
        
        % update Lagrangian gradient
        Gl = g-Jw'*y_new(working)+rho_p*Jw'*(cw-s(working));
        
        % update counter
        k = k+1;
        
        % output iteration parameters
        cv = c(c < -tolerance);
        if k <= 20
            str1 = sprintf (' %3g %4g %7.2f ', k, nf, alpha_k);
            if norm(cv) < tolerance && norm(Gl) < tolerance
                str2 = sprintf ('   (%7.1e)   (%7.1e)    %2u    %7.2f', norm(cv), norm(Gl), l, rho_p);
            elseif norm(cv) < tolerance && norm(Gl) >= tolerance
                str2 = sprintf ('   (%7.1e)    %7.1e     %2u    %7.2f', norm(cv), norm(Gl), l, rho_p);
            elseif norm(cv) >= tolerance && norm(Gl) < tolerance
                str2 = sprintf ('    %7.1e    (%7.1e)    %2u    %7.2f', norm(cv), norm(Gl), l, rho_p);
            elseif norm(cv) >= tolerance && norm(Gl) >= tolerance
                str2 = sprintf ('    %7.1e     %7.1e     %2u    %7.2f', norm(cv), norm(Gl), l, rho_p);
            end
            disp([str1 str2])
        end

    end
    
    % display optimization parameters
    fprintf('\n-> time elapsed: %0.4f\n\n',toc(start_time))
    fprintf('-> number of method iterations: %u\n\n',k)
    fprintf('-> number of function evaluations: %u\n\n',nf)
    fprintf('-> optimizing value of x:\n\n')
    disp(x)
    fprintf('-> optimizing value of y:\n\n')
    disp(y)
    fprintf('-> objective optimum: %0.4f\n\n',f)
    fprintf('-> norm of Lagrangian gradient: %0.4f\n\n',norm(Gl))

end



