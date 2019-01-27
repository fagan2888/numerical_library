% Multivariate (in primal/dual variables) augmented Lagrangian
% SQP merit-function method with modified pos. def. Hessian approx.
% and backtracking line-search.
%
%   1. Solves inequality constrained optimization problem
% 
%           minimize f(x)
%         subject to c(x) >= 0
%
%      with a multivariate augmented Lagrangian serving as a merit-function
%      used for an SQP method.
%
%   2. Uses pos. def. approximation (B) to Lagrangian Hessian (Hl)
%      based on two stage modification:
%
%       stage 1: ensures min. dist. to sing. of reduced Hess. Z'*Hl*Z
%
%       stage 2: ensures pk'*B*pk > 0
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
%        J(acobian matrix of constraint gradients),
%        H(essian of Lagrangian),
%        Hc (Hessian of constraint vector, in general a 3-D matrix)
%       ]
%
%   5. The variable sigma defines the minimum distance from
%      singularity of the reduced Hessian:
%
%           sigma = norm(Z'*H*Z)*beta*eps
%
%      where eps is the machine epsilon and beta is a scaling
%      factor that ensures acceptable error.
%
%   6. The variable 'sol' can be used to select method of KKT solution:
%
%       'w'   : general method -> J[Z Y] = [0 W], Y'Z=0
%       'tq'  : triangular factorization -> J[Z Y] = [0 T], Y'Z=0 (default)
%       'pbs' : basis factorization -> J[Z Y] = [B S], Y'Z=0
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


function [x_new,f_new,g_new] = sqpla_mod_ie(f_name,x,y,sol)

    % define optimization parameters
    nf = 0;                 % init. func. eval. counter
    k = 0;                  % init. method eval. counter
    gamma_c = 1/2;          % line-search contraction factor
    eta_s = 1/4;            % sufficient decrease ratio
    rho_p = 0;              % merit/model-func. penalty parameter
    sigma_b = 0;
    
    % start timer
    start_time = tic;

    % get initial values
    [f,c,g,J,H,~,Hc] = feval(f_name,x);nf=nf+1;
    
    % get system sizes
    n = numel(x);
    m = numel(c);

    % initialize multiplier estimate
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
    beta = 10e03;
    iter_limit = 100;
    
    % construct gradient and Hessian
    Gl = [g-J'*y+rho_p*J'*(c-s);-(c-s)];
    Hl = H - threeD_mult(y,Hc);
    
    % run SQP method w/ modified Hessian routine
    while (norm(Gl) >= tolerance) && (k < iter_limit)

        % output iteration parameters
        cv = c(c < -tolerance);
        if k == 0
            fprintf('\n')
            fprintf(' itn   nf    step    violation    optimal   sigma_b    rho_p\n')
            fprintf('------------------------------------------------------------\n')
            str1 = sprintf (' %3g %4g       - ', k, nf);
            if norm(cv) < tolerance && norm(Gl) < tolerance
                str2 = sprintf ('   (%7.1e)   (%7.1e)    %5.2f  %7.2f', norm(cv), norm(Gl), sigma_b, rho_p);
            elseif norm(cv) < tolerance && norm(Gl) >= tolerance
                str2 = sprintf ('   (%7.1e)    %7.1e   %5.2f    %7.2f', norm(cv), norm(Gl), sigma_b, rho_p);
            elseif norm(cv) >= tolerance && norm(Gl) < tolerance
                str2 = sprintf ('    %7.1e   (%7.1e)     %5.2f  %7.2f', norm(cv), norm(Gl), sigma_b, rho_p);
            elseif norm(cv) >= tolerance && norm(Gl) >= tolerance
                str2 = sprintf ('    %7.1e     %7.1e   %5.2f    %7.2f', norm(cv), norm(Gl), sigma_b, rho_p);
            end
            disp([str1 str2])
        end

        % get e.value matrix of H = LDL' decomposition
        Z = null(J);
        [L,D,P] = ldl(Z'*Hl*Z);
        [V,S] = eig(D);
        diaS = diag(S);

        % get sigma
        sigma = norm(Z'*Hl*Z)*beta*eps;

        % stage 1 modification: dist. to sing. of reduced Hessian
        if min(diaS) < sigma
            
            % modify eigenvalues of H = LDL decomposition
            for i = 1:n-m
                l_bar = abs(diaS(i));
                diaS(i) = max(l_bar,sigma);
            end
            S = diag(diaS);

            % create pos. def. Hessian approximation
            E1 = Z*(P*L*V*S*V'*L'*P - Z'*Hl*Z)*Z';
            B1 = Hl + E1;

        else
            B1 = Hl;
        end

        % solve QP
        pk = quad_prog(J,-c,g,B1,[],0);
        rk = J*pk+c-s;
        
        % stage 2 modification: descent direction
        sigma_b = 0;
        B = B1 + sigma_b*(J'*J);
        if pk'*B*pk <= sigma && norm(pk) > tolerance
            sigma_b = -(pk'*B1*pk)/((rk-c+s)'*(rk-c+s));
            B = B1 + sigma_b*(J'*J);
            while pk'*B*pk <= 0
                sigma_b = 2*sigma_b;
                B = B1 + sigma_b*(J'*J);
            end
        end
        
        % get working sets
        ch = c + J*pk;
        working = ch < tolerance;
        Jw = J(working,:);
        cw = c(working);
        
        % solve KKT system
        [~,yt] = solveKKT(g,B,cw,Jw,sol);
        y_new = zeros(m,1);y_new(working) = yt;
        
        % solve KKT
%         y_new = linsolve(J',g+B*pk);
        qk = y_new - y;
        
        % update penalty parameter
        dphi = (1/2)*(pk'*B*pk)+2*qk'*(c-s)-rk'*y_new;
        if dphi > tolerance
            rho_h = 0;
        else
            rho_h = (rk'*y_new-2*qk'*(c-s))/((c-s+tolerance)'*(c-s+tolerance));
        end
        rho_p = max(rho_h+tolerance,4*rho_p);
        
        % construct merit function
        fp = f-y'*(c-s)+(rho_p/2)*((c-s)'*(c-s));

        % initialize alpha
        alpha_k = 1;
        
        % get initial candidate argument update
        x_new = x + alpha_k*pk;
        s_new = s+alpha_k*rk;
        [f_new,c_new,g_new,J_new,H_new,~,Hc_new] = feval(f_name,x_new);nf=nf+1;

        % construct updated merit function and Lagrangian Hessian
        fp_new = f_new-y_new'*(c_new-s_new)+(rho_p/2)*((c_new-s_new)'*(c_new-s_new));

        % initialize line-search increment
        j = 0;
        
        % perform line-search
        while fp-fp_new < eta_s*alpha_k*(pk'*B*pk-rk'*y_new+(rho_p*(c-s)+2*qk)'*(c-s)) && j < iter_limit/2

            alpha_k = gamma_c*alpha_k;
            x_new = x + alpha_k*pk;
            y_new = y + alpha_k*qk;
            s_new = s+alpha_k*rk;
            [f_new,c_new,g_new,J_new,H_new,~,Hc_new] = feval(f_name,x_new);nf=nf+1;
            fp_new = f_new-y_new'*(c_new-s_new)+(rho_p/2)*((c_new-s_new)'*(c_new-s_new));
            j = j + 1;

        end
        Hl_new = H_new - threeD_mult(y_new,Hc_new);
        
        % update values
        x = x_new;
        y = y_new;
        f = f_new;
        g = g_new;
        c = c_new;
        J = J_new;
        Hl = Hl_new;
        
        % update slack variables
        if rho_p == 0
            s = max(0,c);
        else
            s = max(0,c-y/rho_p);
        end
        
        % update Lagrangian gradient
        Gl = [g-J'*y+rho_p*J'*(c-s);-(c-s)];
        
        % update counters
        k = k+1;
        
        % output iteration parameters
        cv = c(c < -tolerance);
        if k <= 20
            str1 = sprintf (' %3g %4g %7.2f ', k, nf, alpha_k);
            if norm(cv) < tolerance && norm(Gl) < tolerance
                str2 = sprintf ('   (%7.1e)   (%7.1e)  %5.2f    %7.2f', norm(cv), norm(Gl), sigma_b, rho_p);
            elseif norm(cv) < tolerance && norm(Gl) >= tolerance
                str2 = sprintf ('   (%7.1e)    %7.1e   %5.2f    %7.2f', norm(cv), norm(Gl), sigma_b, rho_p);
            elseif norm(cv) >= tolerance && norm(Gl) < tolerance
                str2 = sprintf ('    %7.1e   (%7.1e)   %5.2f    %7.2f', norm(cv), norm(Gl), sigma_b, rho_p);
            elseif norm(cv) >= tolerance && norm(Gl) >= tolerance
                str2 = sprintf ('    %7.1e     %7.1e   %5.2f    %7.2f', norm(cv), norm(Gl), sigma_b, rho_p);
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

