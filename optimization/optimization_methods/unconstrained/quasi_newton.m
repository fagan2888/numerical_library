% Quasi-Newton method with Armijo line-search and BFGS update.
%
%   1. Implements a Hessian-free quasi-Newton method, using
%      the BFGS algorithm to update the forward-difference Hessian (FDH).
%
%   2. Uses pos. def. approximation (B) to FDH (H)
%      based on modifying the eigenvalues of the H = LDL' decomp.
%      performed by a symmetric indefinite solver.
%
%   3. Performs Armijo line-search for acceptable objective reduction.
%
%   4. Calls objective function, which must return
%
%           [f(unction),g(radient),H(essian)]
%
%      (the Hessian is unused in this method).
%
%   5. The variable sigma defines the minimum distance from
%      singularity of the pos. def. approximation:
%
%           sigma = norm(H)*beta*eps
%
%      where eps is the machine epsilon and beta is a scaling
%      factor that ensures acceptable error.

% Jeremy Orosco
% MATH 271A, Fall 2015
% Coimbra Energy Group
% Department of Mechanical
% and Aerospace Engineering
% University of California, San Diego

function [x_new,f_new,g_new] = quasi_newton(f_name,x)

    % get system size
    n = numel(x);

    % define optimization parameters
    k = 0;
    nf = 0;
    gamma_c = 1/2;
    eta_s = 1/4;
    
    % define limits
    tolerance = 10*sqrt(eps);
    beta = 10e03;
    iter_limit = 10e03;
    
    % start timer
    start_time = tic;

    % get initial values
    [f,g,~] = feval(f_name,x);
    
    % initialize Hessian with identity
    H = eye(n);
    
    % run quasi-Newton method
    while norm(g) >= tolerance && k < iter_limit
    
        % get e.value matrix of H = LDL' decomposition
        [L,D,P] = ldl(H);
        [V,S] = eig(D);
        diaS = diag(S);

        % get sigma
        sigma = norm(H)*beta*eps;
        
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
            
            B = H;
            
        end

        % get descent direction
        p = -linsolve(B,g);

        % initialize alpha
        alpha_k = 1;

        % get initial candidate argument update
        x_new = x + alpha_k*p;
        [f_new,g_new,~] = feval(f_name,x_new);

        % perform Armijo line search
        while f_new > f + alpha_k*eta_s*g'*p

            alpha_k = gamma_c*alpha_k;
            x_new = x + alpha_k*p;
            [f_new,g_new,~] = feval(f_name,x_new);
            j = j + 1;
            nf = nf + 1;

        end

        % compute FDH update
        dk = x_new - x;
        yk = g_new - g;
        H_new = B+(yk*yk')/(yk'*dk)-(B*(dk*dk')*B')/(dk'*B*dk); % BFGS

        % update values
        x = x_new;
        f = f_new;
        g = g_new;
        H = H_new;
        
        % update counters
        k = k + 1;
        nf = nf + 1;
    
    end
    
    % display elapsed time
    disp(toc(start_time))
    
    % display optimization parameters
    space = sprintf('\n');
    fprintf('\nnumber of iterations:\n')
    disp(space)
    disp(k)
    fprintf('\nminimizing value of x:\n')
    disp(space)
    disp(x)
    fprintf('\nfunction minimum:\n')
    disp(space)
    disp(f)
    fprintf('\ngradient at minimizing value of x:\n')
    disp(space)
    disp(g)

end

