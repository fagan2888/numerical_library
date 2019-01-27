% Multivariate (in primal/dual variables) augmented Lagrangian
% SQP merit-function method with modified pos. def. Hessian approx.
% and backtracking line-search.
%
%   1. Solves constrained optimization with multivariate augmented
%      Lagrangian serving as a merit-function for an SQP method.
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


function [x_new,f_new,g_new] = sqpla_mod_e(f_name,x,y,sol)

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
    if (nargin < 4)
        if (nargin < 3)
            y = 10e-03*ones(m,1);
        end
        sol = 'tq';
    end
    
    % define limits
    tolerance = 10e02*sqrt(eps);
    beta = 10e03;
    iter_limit = 100;
    
    % construct gradient and Hessian
    Gl = [g-J'*y;-c];
    Hl = H - threeD_mult(y,Hc);
    
    % run SQP method w/ modified Hessian routine
    while (norm(Gl) >= tolerance) && (k < iter_limit)

        % output iteration parameters
        if k == 0
            fprintf('\n')
            fprintf(' itn   nf    step    feasible     optimal   sigma_b    rho_p\n')
            fprintf('------------------------------------------------------------\n')
            str1 = sprintf (' %3g %4g       - ', k, nf);
            if norm(c) < tolerance && norm(Gl) < tolerance
                str2 = sprintf ('   (%7.1e)   (%7.1e)    %5.2f  %7.2f', norm(c), norm(Gl), sigma_b, rho_p);
            elseif norm(c) < tolerance && norm(Gl) >= tolerance
                str2 = sprintf ('   (%7.1e)    %7.1e   %5.2f    %7.2f', norm(c), norm(Gl), sigma_b, rho_p);
            elseif norm(c) >= tolerance && norm(Gl) < tolerance
                str2 = sprintf ('    %7.1e   (%7.1e)     %5.2f  %7.2f', norm(c), norm(Gl), sigma_b, rho_p);
            elseif norm(c) >= tolerance && norm(Gl) >= tolerance
                str2 = sprintf ('    %7.1e     %7.1e   %5.2f    %7.2f', norm(c), norm(Gl), sigma_b, rho_p);
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

        % solve KKT system with preferred factorization
        [pk,y_new] = solveKKT(g,B1,c,J,sol);
        
        % stage 2 modification: descent direction
        sigma_b = 0;
        B = B1 + sigma_b*(J'*J);
        if pk'*B*pk <= sigma
            sigma_b = -(pk'*B1*pk)/(c'*c);
            B = B1 + sigma_b*(J'*J);
            while pk'*B*pk <= 0
                sigma_b = 2*sigma_b;
                B = B1 + sigma_b*(J'*J);
            end
        end
        % correct singular sigma_b for real return
        % IS THIS VALID?
        if abs(sigma_b) > 10e6
            sigma_b = sign(sigma_b)*10e6;
        end
        y_new = y_new - sigma_b*c;
        qk = y_new - y;
        
        % update penalty parameter
        dphi = -((1/2)*(pk'*B*pk)+2*qk'*c+rho_p*(c'*c));
        if dphi > 0
            if (1/2)*pk'*B*pk <= -2*qk'*c
                rho_h = 2*(norm(qk)/norm(c));
            else
                rho_h = 0;
            end
            rho_p = max(rho_h+tolerance,2*rho_p);
        end
        
        % construct merit function
        fp = f - y'*c + (rho_p/2)*(c'*c);

        % initialize alpha
        alpha_k = 1;
        
        % get initial candidate argument update
        x_new = x + alpha_k*pk;
        [f_new,c_new,g_new,J_new,H_new,~,Hc_new] = feval(f_name,x_new);nf=nf+1;

        % construct updated merit function and Lagrangian Hessian
        fp_new = f_new - y_new'*c_new + (rho_p/2)*(c_new'*c_new);

        % initialize line-search increment
        j = 0;
        
        % perform line-search
        while fp - fp_new < eta_s*alpha_k*(pk'*B*pk+2*qk'*c+rho_p*(c'*c)) && j < iter_limit/2

            alpha_k = gamma_c*alpha_k;
            x_new = x + alpha_k*pk;
            y_new = y + alpha_k*qk;
            [f_new,c_new,g_new,J_new,H_new,~,Hc_new] = feval(f_name,x_new);nf=nf+1;
            fp_new = f_new - y_new'*c_new + (rho_p/2)*(c_new'*c_new);
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
        Gl = [g-J'*y;-c];

        % update counters
        k = k+1;

        % output iteration parameters
        if k <= 20
            str1 = sprintf (' %3g %4g %7.2f ', k, nf, alpha_k);
            if norm(c) < tolerance && norm(Gl) < tolerance
                str2 = sprintf ('   (%7.1e)   (%7.1e)  %5.2f    %7.2f', norm(c), norm(Gl), sigma_b, rho_p);
            elseif norm(c) < tolerance && norm(Gl) >= tolerance
                str2 = sprintf ('   (%7.1e)    %7.1e   %5.2f    %7.2f', norm(c), norm(Gl), sigma_b, rho_p);
            elseif norm(c) >= tolerance && norm(Gl) < tolerance
                str2 = sprintf ('    %7.1e   (%7.1e)   %5.2f    %7.2f', norm(c), norm(Gl), sigma_b, rho_p);
            elseif norm(c) >= tolerance && norm(Gl) >= tolerance
                str2 = sprintf ('    %7.1e     %7.1e   %5.2f    %7.2f', norm(c), norm(Gl), sigma_b, rho_p);
            end
            disp([str1 str2])
        end

    end
    
    % temporary output format adjustment
    format long
    
    % display optimization parameters
    space = sprintf('\n');
    fprintf('\ntime elapsed:')
    disp(toc(start_time))
    fprintf('number of method iterations:')
    disp(k)
    fprintf('number of function evaluations:')
    disp(nf)
    fprintf('optimizing value of x:')
    disp(space)
    disp(x)
    fprintf('optimizing value of y:')
    disp(space)
    disp(y)
    fprintf('objective optimum:')
    disp(f_new)
    fprintf('Lagrangian gradient at optimum:')
    disp(space)
    disp(Gl)
    
    % return to original formatting
    format short

end

