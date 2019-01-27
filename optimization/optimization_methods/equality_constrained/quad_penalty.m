% Penalty function with Modified Newton's method and Armijo line search.
%
%   1. Solves constrained optimization with penalty function method.
%
%   2. Uses pos. def. approximation (B) to Hessian (H)
%      based on modifying the eigenvalues of the H = LDL' decomp.
%      performed by a symmetric indefinite solver.
%
%   3. Performs Armijo line-search for acceptable objective reduction.
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
%      singularity of the pos. def. approximation:
%
%           sigma = norm(H)*beta*eps
%
%      where eps is the machine epsilon and beta is a scaling
%      factor that ensures acceptable error.
%
%   6. The input parameter opt with possible values 'max' and 'min'
%      indicates the type of optima to be located. The default value is
%      'min'.
%
%   7. Uses gradient of penalty function to determine convergence of
%      penalty function (inner) loop. Uses absolute change in objective
%      function to determine convergence of method (outer) loop.

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


function [x_new,f_new,g_new] = quad_penalty(f_name,x)

    % get system size
    n = numel(x);

    % define optimization parameters
    nf = 0; % function evaluation counter
    k = 0;
    l = 0; % method iteration counter
    gamma_c = 1/2;
    eta_s = 1/4;
    gamma = 10; % expansion factor for rho
    rho = 1; % initial value
    
    % define limits
    tolerance = 10e01*sqrt(eps);
    beta = 10e03;
    iter_limit = 100;
    
    % start timer
    start_time = tic;
    
    % initialize df
    df = 1; % method convergence parameter
    
    % run penalty function loop
    while (df >= tolerance) && (l < iter_limit)

        % get initial values
        [f,c,g,J,H,~,Hc] = feval(f_name,x);nf=nf+1;
        f_new = f;
        
        % construct penalty function
        fp = f + (rho/2)*(c'*c);
        gp = g + rho*(J'*c);
        Hp = H + rho*threeD_mult(c,Hc) + rho*(J'*J);
    
        % run modified Newton's method
        while (norm(gp) >= tolerance) && (k < iter_limit)

            % get e.value matrix of H = LDL' decomposition
            [L,D,P] = ldl(Hp);
            [V,S] = eig(D);
            diaS = diag(S);

            % get sigma
            sigma = norm(Hp)*beta*eps;

            % initialize line-search increment
            j = 0;

            % check for min(eig. val.) > sigma, modify if necessary
            if min(diaS) < sigma

                % modify eigenvalues of H = LDL decomposition
                for i = 1:n
                    l_bar = abs(diaS(i));
                    diaS(i) = max(l_bar,sigma);
                end
                S = diag(diaS);

                % create pos. def. Hessian approximation
                B = P*L*V*S*V'*L'*P';

            else

                B = Hp;

            end

            % get descent direction
            pk = linsolve(B,-gp);

            % initialize alpha
            alpha_k = 1;

            % get initial candidate argument update
            x_new = x + alpha_k*pk;
            [f_new,c,g_new,J,H,~,Hc] = feval(f_name,x_new);nf=nf+1;

            % construct penalty function
            fp_new = f_new + (rho/2)*(c'*c);
            gp_new = g_new + rho*(J'*c);
            Hp_new = H + rho*threeD_mult(c,Hc) + rho*(J'*J);

            % perform Armijo line search
            while (fp_new > fp + alpha_k*eta_s*(gp'*pk)) && (j < iter_limit/2)

                alpha_k = gamma_c*alpha_k;
                x_new = x + alpha_k*pk;
                [f_new,c,g_new,J,H,~,Hc] = feval(f_name,x_new);nf=nf+1;
                fp_new = f_new + (rho/2)*(c'*c);
                gp_new = g_new + rho*(J'*c);
                Hp_new = H + rho*threeD_mult(c,Hc) + rho*(J'*J);
                j = j + 1;

            end

            % update values
            x = x_new;
            fp = fp_new;
            gp = gp_new;
            Hp = Hp_new;

            % update counters
            k = k+1;

        end
        
        % increase penalty parameter
        rho = gamma*rho;
        
        % calculate change in objective
        df = abs(f-f_new);
        
        % update counter
        l = l + 1;
    
    end
    
    % display elapsed time
    disp(toc(start_time))
    
    % display optimization parameters
    space = sprintf('\n');
    fprintf('\nnumber of method iterations:\n')
    disp(space)
    disp(l)
    fprintf('\nnumber of function evaluations:\n')
    disp(space)
    disp(nf)
    fprintf('\noptimizing value of x:\n')
    disp(space)
    disp(x)
    fprintf('\nobjective optimum:\n')
    disp(space)
    disp(f_new)
    fprintf('\npenalty function gradient at optimum:\n')
    disp(space)
    disp(gp)

end


% vector / 3D matrix multiplication
function product = threeD_mult(v,M)

    n = size(M,1);
    m = size(M,3);

    sum = zeros(n);
    for i = 1:m
        sum = sum + v(i)*M(:,:,i);
    end
    
    product = sum;

end

