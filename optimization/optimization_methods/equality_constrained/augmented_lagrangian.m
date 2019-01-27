% Augmented Lagrangian method using second-order Lagrange multiplier
% update and with a Modified Newton's method using an Armijo line-search.
%
%   1. Solves constrained optimization with augmented Lagrangian method.
%
%   2. Updates to multiplier estimates are of second-order, so that
%      overall method retains Q-quadratic convergence.
%
%   3. Uses pos. def. approximation (B) to Hessian (H)
%      based on modifying the eigenvalues of the H = LDL' decomp.
%      performed by a symmetric indefinite solver.
%
%   4. Performs Armijo line-search for acceptable objective reduction.
%
%   5. Calls objective function, which must return
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
%   6. The variable sigma defines the minimum distance from
%      singularity of the pos. def. approximation:
%
%           sigma = norm(H)*beta*eps
%
%      where eps is the machine epsilon and beta is a scaling
%      factor that ensures acceptable error.
%
%   7. The input parameter opt with possible values 'max' and 'min'
%      indicates the type of optima to be located. The default value
%      is 'min'.
%
%   8. Uses gradient of augmented Lagrangian with current multiplier
%      estimate to determine convergence of augmented Lagrangian (inner)
%      loop. Uses gradient of augmented Lagrangian with multiplier update
%      to determine convergence of method (outer) loop.
%
%   9. A bound is enforced on the Lagrange multiplier estimates, so that
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


function [x,f,g] = augmented_lagrangian(f_name,x,y,opt)

    % set optimum type
    if (nargin < 4) || strcmp(opt,'min')
        if (nargin < 3)
            y = 0;
        end
        opt = 1;
    elseif strcmp(opt,'max')
        opt = -1;
    end

    % get system size
    n = numel(x);
    m = numel(y);

    % define optimization parameters
    nf = 0;             % initialize function evaluation counter
    k = 0;              % initialze aug. Lagr. funct. min. counter
    l = 0;              % initialize method iteration counter
    gamma_c = 1/2;      % line-search contraction factor
    eta_s = 1/4;        % sufficient decrease ratio
    gamma = 10;         % penalty parameter multiplier
    rho = 1;            % -> inf. penalty paramter
    e = ones(m,1);      % unit vector
    
    % define limits
    tolerance = 10e02*sqrt(eps);    % convergence
    beta = 10e03;                   % singularity distance factor
    iter_limit = 100;               % loop-break limit
    y_max = 10^5;                   % Lagrange mult. estimate bound
    
    % start timer
    start_time = tic;
    
    % get initial values
    [f,c,g,J,H,~,Hc] = feval(f_name,x);nf=nf+1;
    
    % construct augmented Lagrangian
    piy = y - rho*c;
    La = opt*f - c'*(y - (rho/2)*c);
    Ga = opt*g - J'*piy;
    Hp = opt*H - threeD_mult(piy,Hc);
    Ha = Hp + rho*(J'*J);
    
    % run method loop
    while(norm(Ga) >= tolerance) && (l < iter_limit)
    
        % run modified Newton's method
        while (norm(Ga) >= tolerance) && (k < iter_limit)

            % get e.value matrix of Kl = LDL' decomposition
            [L,D,P] = ldl(Ha);
            [V,S] = eig(D);
            diaS = diag(S);

            % get sigma
            sigma = norm(Ha)*beta*eps;

            % check for min(eig. val.) > sigma, modify if necessary
            if min(diaS) < sigma

                % modify eigenvalues of H = LDL decomposition
                for i = 1:n
                    l_bar = abs(diaS(i));
                    diaS(i) = max(l_bar,sigma);
                end
                S = diag(diaS);

                % create pos. def. Hessian approximation
                Ba = P*L*V*S*V'*L'*P';

            else

                Ba = Ha;

            end

            % get descent direction
            pk = linsolve(Ba,-Ga);

            % initialize alpha
            alpha_k = 1;

            % get initial candidate argument update
            x_new = x + alpha_k*pk;
            [f,c,g,J,H,~,Hc] = feval(f_name,x_new);nf=nf+1;

            % construct aug. Lagrangian function
            piy = y - rho*c;
            La_new = opt*f - c'*(y - (rho/2)*c);
            Ga_new = opt*g - J'*piy;
            Hp_new = opt*H - threeD_mult(piy,Hc);
            Ha_new = Hp_new + rho*(J'*J);

            % initialize line-search counter
            j = 0;

            % perform Armijo line search
            while (La_new > La + alpha_k*eta_s*(Ga'*pk)) && (j < iter_limit/2)

                alpha_k = gamma_c*alpha_k;
                x_new = x + alpha_k*pk;
                [f,c,g,J,H,~,Hc] = feval(f_name,x_new);nf=nf+1;
                piy = y - rho*c;
                La_new = opt*f - c'*(y - (rho/2)*c);
                Ga_new = opt*g - J'*piy;
                Hp_new = opt*H - threeD_mult(piy,Hc);
                Ha_new = Hp_new + rho*(J'*J);
                j = j + 1;

            end

            % update values
            x = x_new;
            La = La_new;
            Ga = Ga_new;
            Ha = Ha_new;
            
            % update counter
            k = k+1;

        end
        
        % increase penalty parameter
        rho = gamma*rho;

        % update Lagrange multiplier
        py = linsolve((J*(Ha^-1)*J'),-c);   % second-order update
        y = y + py;                         % get estimate update
        y = max(-y_max*e,min(y,y_max*e));   % bound estimate

        % construct augmented Lagrangian
        piy = y - rho*c;
        La = opt*f - c'*(y - (rho/2)*c);
        Ga = opt*g - J'*piy;
        Hp = opt*H - threeD_mult(piy,Hc);
        Ha = Hp + rho*(J'*J);
        
        % update counter
        l = l + 1;

    end
    
    % display optimization parameters
    space = sprintf('\n');
    fprintf('\ntime taken:\n')
    disp(space)
    disp(toc(start_time))
    fprintf('\nnumber of method iterations:\n')
    disp(space)
    disp(l)
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

