% L1-SQP merit-function method with quasi-Newton method using BFGS update
% and backtracking line-search. Includes problem scale based method
% preference for efficient solution of the corresponding KKT equations.
%
%   1. Solves constrained optimization with L1 penalty function serving
%      as a merit-function for an SQP method.
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
% MATH 271B, Fall 2016
% Coimbra Energy Group
% Department of Mechanical
% and Aerospace Engineering
% University of California, San Diego

% Copyright 2016 Jeremy Orosco
% This work is licensed under the
% Creative Commons Attribution-NonCommercial-ShareAlike 4.0
% International License. To view a copy of this license,
% visit http://creativecommons.org/licenses/by-nc-sa/4.0/.


function [x_new,f_new,g_new] = sl1qp_quasi(f_name,x,y,sol)

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

    % initialize multiplier estimate
    if (nargin < 4)
        if (nargin < 3)
            y = 10e-03*ones(m,1);
        end
        sol = 'tq';
    end
    
    % define limits
    tolerance = 10e02*sqrt(eps);
    iter_limit = 100;
    
    % construct gradient and Hessian
    Gl = g - J'*y;
    B = eye(n);
    
    % run L1-SQP method w/ modified Hessian routine
    while (norm(Gl) >= tolerance) && (k < iter_limit)

        % output iteration parameters
        if k == 0
            fprintf('\n')
            fprintf(' itn   nf    step    feasible     optimal    nz     rho\n')
            fprintf('---------------------------------------------------------\n')
            str1 = sprintf (' %3g %4g       - ', k, nf);
            if norm(c) < tolerance && norm(Gl) < tolerance
                str2 = sprintf ('   (%7.1e)   (%7.1e)    %u', norm(c), norm(Gl), l);
            elseif norm(c) < tolerance && norm(Gl) >= tolerance
                str2 = sprintf ('   (%7.1e)    %7.1e     %u', norm(c), norm(Gl), l);
            elseif norm(c) >= tolerance && norm(Gl) < tolerance
                str2 = sprintf ('    %7.1e   (%7.1e)     %u', norm(c), norm(Gl), l);
            elseif norm(c) >= tolerance && norm(Gl) >= tolerance
                str2 = sprintf ('    %7.1e     %7.1e     %u', norm(c), norm(Gl), l);
            end
            disp([str1 str2])
        end

        % solve KKT system with preferred factorization
        [pk,y,pn] = solveKKT(g,B,c,J,sol);

        % update merit/model function penalty parameter
        rho_p = max(norm(y,inf)+tolerance,rho_p);

        % construct merit function
        fp = f + (rho_p/2)*norm(c,1);

        % initialize alpha
        alpha_k = 1;
        
        % get initial candidate argument update
        x_new = x + alpha_k*pk;
        [f_new,c_new,g_new,J_new,~,~,~] = feval(f_name,x_new);nf=nf+1;
        sk = -J'*((J*J')^-1)*c_new;

        % construct updated merit function and Lagrangian Hessian
        fp_new = f_new + (rho_p/2)*norm(c_new,1);
        Gl_new = g_new - J_new'*y;

        % initialize line-search increment
        j = 0;
        
        % perform Maratos safeguarded line-search
        while fp - fp_new < eta_s*alpha_k*(rho_p*norm(c,1)-g'*pk) && j < iter_limit/2

            if j == 0
                pk_new = pk + sk;
                x_new = x + pk_new;
                [f_new,c_new,g_new,J_new,~,~,~] = feval(f_name,x_new);nf=nf+1;
                fp_new = f_new + (rho_p/2)*norm(c_new,1);
                j = j + 1;
                if fp - fp_new < eta_s*alpha_k*(rho_p*norm(c,1)-g'*pk_new)
                    break
                end
            else
                alpha_k = gamma_c*alpha_k;
                pk_new = alpha_k*pk + sk;
                x_new = x + pk_new;
                [f_new,c_new,g_new,J_new,~,~,~] = feval(f_name,x_new);nf=nf+1;
                fp_new = f_new + (rho_p/2)*norm(c_new,1);
                j = j + 1;
                if fp - fp_new < eta_s*alpha_k*(rho_p*norm(c,1)-g'*pk_new)
                    break
                end
            end

        end
        
        % compute FDH update
        dk = x_new - x;
        yk = Gl_new - Gl;

        % compute BFGS update
        sigma = alpha_k*(1-eta_h)*pk'*B*pk;
        if yk'*dk >= sigma
            B = B+(yk*yk')/(yk'*dk)-(B*(dk*dk')*B')/(dk'*B*dk);
        else
            z = x + alpha_k*pn;
            dk = x_new - z;
            [~,~,g_z,J_z,~,~,~] = feval(f_name,z);nf=nf+1;
            Gl_z = g_z - J_z'*y;
            yk = Gl_new - Gl_z;
            if yk'*dk >= sigma
                B = B+(yk*yk')/(yk'*dk)-(B*(dk*dk')*B')/(dk'*B*dk);
            end
            l = l + 1;
        end

        % update values
        x = x_new;
        f = f_new;
        c = c_new;
        g = g_new;
        J = J_new;
        Gl = Gl_new;

        % update counters
        k = k+1;

        % output iteration parameters
        if k <= 20
            str1 = sprintf (' %3g %4g %7.2f ', k, nf, alpha_k);
            if norm(c) < tolerance && norm(Gl) < tolerance
                str2 = sprintf ('   (%7.1e)   (%7.1e)    %u    %.2f', norm(c), norm(Gl), l, rho_p);
            elseif norm(c) < tolerance && norm(Gl) >= tolerance
                str2 = sprintf ('   (%7.1e)    %7.1e     %u    %.2f', norm(c), norm(Gl), l, rho_p);
            elseif norm(c) >= tolerance && norm(Gl) < tolerance
                str2 = sprintf ('    %7.1e   (%7.1e)     %u    %.2f', norm(c), norm(Gl), l, rho_p);
            elseif norm(c) >= tolerance && norm(Gl) >= tolerance
                str2 = sprintf ('    %7.1e     %7.1e     %u    %.2f', norm(c), norm(Gl), l, rho_p);
            end
            disp([str1 str2])
        end

    end
    
    % display optimization parameters
    space = sprintf('\n');
    fprintf('\ntime elapsed:\n')
    disp(space)
    disp(toc(start_time))
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
    disp(Gl)

end


% solve KKT system efficiently
function [pk,y,pn] = solveKKT(g,B,c,J,opt)

    if strcmp(opt,'tq') == 1
        
        % get system size
        n = numel(g);
        m = numel(c);

        % construct change of variables (JQ = [0 T])
        [Qt,Tt] = qr(J'); % get transpose variables
        P = fliplr(eye(n)); % create permutation matrix
        T = (P*Tt)';T = T(:,n-m+1:end); % get T factor
        Q = Qt*P; % get Q factor
        Z = Q(:,1:n-m); % nullspace of J
        Y = Q(:,n-m+1:end); % rangespace of J'

        % construct QP subproblem solution
        pn = Y*linsolve(T,-c); % feasibility step
        pt = Z*linsolve(Z'*B*Z,-Z'*(g+B*pn)); % minimization step
        pk = pt + pn; % overall solution step

        % construct optimal QP subproblem multiplier
        y = linsolve(T',Y'*(g+B*pk));
        
    elseif strcmp(opt,'pbs') == 1
        
        % get system size
        n = numel(g);
        m = numel(c);
        
        % construct change of variables (JP = [Bj S])
        [~,~,Pt] = lu(J');P = Pt';
        BjS = J*P;
        Bj = BjS(:,1:m);
        S = BjS(:,m+1:end);
        Z = P*[-linsolve(Bj,S);eye(n-m)];
        Y = P*[eye(m);zeros(n-m,m)];
        
        % construct QP subproblem solution
        pn = Y*linsolve(Bj,-c); % feasibility step
        pt = Z*linsolve(Z'*B*Z,-Z'*(g+B*pn)); % minimization step
        pk = pt + pn; % overall solution step
        
        % construct optimal QP subproblem multiplier
        y = linsolve(Bj',Y'*(g+B*pk));
        
    elseif strcmp(opt,'w') == 1
        
        % construct change of variables (JQ = [0 W])
        Z = null(J);
        Y = orth(J');
        W = J*Y;
        
        % construct QP subproblem solution
        pn = Y*linsolve(W,-c); % feasibility step
        pt = Z*linsolve(Z'*B*Z,-Z'*(g+B*pn)); % minimization step
        pk = pt + pn; % overall solution step
        
        % construct optimal QP subproblem multiplier
        y = linsolve(W',Y'*(g+B*pk));
        
    end


end

