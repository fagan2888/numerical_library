% Augmented Lagrangian method using a KKT formulation with a inertia
% conserving LDL' factorization for obtaining a modified KKT matrix.
%
%   1. Solves constrained optimization with augmented Lagrangian method.
%
%   2. Inertia conserving LDL' factorization ensures a solution to the
%      KKT system that provides a descent direction.
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
%   5. The input parameter opt with possible values 'max' and 'min'
%      indicates the type of optima to be located. The default value
%      is 'min'.
%
%   7. Uses gradient of augmented Lagrangian to determine convergence.
%
%   8. A bound is enforced on the Lagrange multiplier estimates, so that
%      for large multiplier estimates the algorithm behaves like a
%      classical penalty method, with a minimizer obtained for rho -> inf.

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


function [x,f,g] = augmented_lagrangian_KKT(f_name,x,ye)

    % define optimization parameters
    nf = 0;             % initialize function evaluation counter
    k = 0;              % initialze method iteration counter
    gamma_c = 1/2;      % line-search contraction factor
    eta_s = 1/4;        % sufficient decrease ratio
    gamma = 10;         % penalty parameter multiplier
    rho = 1;            % -> inf. penalty paramter
    
    % start timer
    start_time = tic;
    
    % get initial values
    [f,c,g,J,H,~,Hc] = feval(f_name,x);nf=nf+1;

    % get system size
    n = numel(x);
    m = numel(c);

    % set optimum type
    if (nargin <3)
        ye = zeros(m,1);
    end
    
    % define limits
    tolerance = 10e02*sqrt(eps);    % convergence
    iter_limit = 100;             % loop-break limit
    y_max = 10^5;                   % Lagrange mult. estimate bound
    e = ones(m,1);                  % unit vector
    
    % construct augmented Lagrangian
    piy = ye - rho*c;
    La = f - c'*(ye - (rho/2)*c);
    Ga = g - J'*piy;
    Hp = H - threeD_mult(piy,Hc);
    
    % define multiplier update criteria
    eps_y = 3/4;                % acceptable approx. mult. convergence
    beta_y = (1/4)*norm(c);     % initialize allowable dist. to constraint
    
    % run modified Newton's method
    while (norm(Ga) >= tolerance) && (k < iter_limit)

        % get eigenvalue modification
        [rho,E_k] = mod_ldl(Hp,J,rho);
        
        % construct KKT system
        Kl = [Hp+E_k, -J';
              J     , (1/rho)*eye(m)];
        Kr = [g;
              c-(1/rho)*ye];

        % get descent direction and estimate update
        mk = linsolve(Kl,-Kr);
        pk = mk(1:n);
        y = mk(n+1:n+m);
        
        % initialize alpha
        alpha_k = 1;

        % get initial candidate argument update
        x_new = x + alpha_k*pk;
        [f,c,g,J,H,~,Hc] = feval(f_name,x_new);nf=nf+1;

        % construct aug. Lagrangian function
        piy = ye - rho*c;
        La_new = f - c'*(ye - (rho/2)*c);
        Ga_new = g - J'*piy;

        % initialize line-search counter
        j = 0;

        % perform Armijo line search
        while La_new > La + alpha_k*eta_s*(Ga'*pk) && j < iter_limit/2

            alpha_k = gamma_c*alpha_k;
            x_new = x + alpha_k*pk;
            [f,c,g,J,H,~,Hc] = feval(f_name,x_new);nf=nf+1;
            piy = ye - rho*c;
            La_new = f - c'*(ye - (rho/2)*c);
            Ga_new = g - J'*piy;
            j = j + 1;

        end

        % update values
        x = x_new;
        Ga = Ga_new;
        
        % update multiplier estimate according to criteria
        if norm(Ga) <= eps_y
           if norm(c) >= beta_y
               rho = gamma*rho;
           end
           ye = max(-y_max*e,min(y,y_max*e));
           beta_y = (1/4)*norm(c);
           eps_y = (1/2)*eps_y;
        end
        
        % reconstruct aug. Lagrangian for next iteration
        piy = ye - rho*c;
        La = f - c'*(ye - (rho/2)*c);
        Ga = g - J'*piy;
        Hp = H - threeD_mult(piy,Hc);
        
        % update counter
        k = k+1;

        % output iteration parameters
        if k <= 20
            if k == 1
                fprintf(' itn   nf    step      feasible    optimal\n')
            end
            str1 = sprintf (' %3g %4g %9.2e ', k, nf, alpha_k);
            if norm(c) < tolerance && norm(Ga) < tolerance
                str2 = sprintf ('   (%7.1e)   (%7.1e)', norm(c), norm(Ga));
            elseif norm(c) < tolerance && norm(Ga) >= tolerance
                str2 = sprintf ('   (%7.1e)    %7.1e', norm(c), norm(Ga));
            elseif norm(c) >= tolerance && norm(Ga) < tolerance
                str2 = sprintf ('    %7.1e   (%7.1e)', norm(c), norm(Ga));
            elseif norm(c) >= tolerance && norm(Ga) >= tolerance
                str2 = sprintf ('    %7.1e     %7.1e', norm(c), norm(Ga));
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
    disp(ye)
    fprintf('\nobjective optimum:\n')
    disp(space)
    disp(f)
    fprintf('\naugmented Lagrangian gradient at optimum:\n')
    disp(space)
    disp(Ga)

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


% inertia conserving symmetric indefinite LDL decomposition
function [rho,E_k] = mod_ldl(H,J,rho)

    n = size(H,1);
    m = size(J,1);
    mu = 1/rho;

    sigma = 1;
    E_k = zeros(n);

    K = [H -J';
         J  mu*eye(m)];

    [~,D,~] = ldl(K);
    [~,S] = eig(D);
    diaS = diag(S);
    
    [peig,neig,zeig] = deal(0);
    for i = 1:m+n
        if diaS(i) > 0
            peig = peig + 1;
        elseif diaS(i) < 0
            neig = neig + 1;
        else
            zeig = zeig +1;
        end
    end
 
    while neig > m || zeig > 0
        
        E_k = sigma*eye(n);

        K = [H+E_k -J';
             J      mu*eye(m)];

        [~,D,~] = ldl(K);
        [~,S] = eig(D);
        diaS = diag(S);
    
        [peig,neig,zeig] = deal(0);
        for i = 1:m+n
            if diaS(i) > 0
                peig = peig + 1;
            elseif diaS(i) < 0
                neig = neig + 1;
            else
                zeig = zeig +1;
            end
        end
        
        sigma = 10*sigma;
        
    end
    
    if norm(E_k) > 0
        rho = 10*rho;
    end

end

