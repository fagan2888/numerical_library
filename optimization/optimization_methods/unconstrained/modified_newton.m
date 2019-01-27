% Modified Newton's method with Armijo line search.
%
%   1. Uses pos. def. approximation (B) to Hessian (H)
%      based on modifying the eigenvalues of the H = LDL' decomp.
%      performed by a symmetric indefinite solver.
%
%   2. Performs Armijo line-search for acceptable objective reduction.
%
%   3. Calls objective function, which must return
%
%           [f(unction),g(radient),H(essian)]
%
%   4. The variable sigma defines the minimum distance from
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

function [x_new,f_new,g_new] = modified_newton(f_name,x)

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

    % get initial values
    [f,g,H] = feval(f_name,x);
    
    % start timer
    start_time = tic;
    
    % run modified Newton's method
    while norm(g) >= tolerance && k < iter_limit
    
        % get e.value matrix of H = LDL' decomposition
        [L,D,P] = ldl(H);
        [V,S] = eig(D);
        diaS = diag(S);

        % get sigma
        sigma = norm(H)*beta*eps;
        
        % initialize line search increment
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
        p = -(B^-1)*g;

        % initialize alpha
        alpha_k = 1;
        
        % get initial candidate argument update
        x_new = x + alpha_k*p;
        [f_new,g_new,H_new] = feval(f_name,x_new);

        % perform Armijo line search
        while f_new > f + alpha_k*eta_s*g'*p

            alpha_k = gamma_c*alpha_k;
            x_new = x + alpha_k*p;
            [f_new,g_new,H_new] = feval(f_name,x_new);
            j = j + 1;
            nf = nf + 1;

        end

        % update values
        x = x_new;
        f = f_new;
        g = g_new;
        H = H_new;

        % update counters
        k = k+1;
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

