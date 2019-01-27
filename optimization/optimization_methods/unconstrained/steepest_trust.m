function [x_new,f_new,g_new] = steepest_trust(f_name,x)

    k = 0;
    delta_k = 1;
    gamma_c = 1/2;
    gamma_e = 2;
    eta_s = 1/10;
    eta_e = 1/4;
    tol = 10*sqrt(eps);
    iter_limit = 10e03;
    
    start_time = tic;

    [f,g,~] = feval(f_name,x);

    while norm(g) > tol && k < iter_limit
        
        sigma = (1/delta_k)*norm(g);
        d_k = -(1/sigma)*g;
        
        x_new = x + d_k;
        
        gTd = g'*d_k;
        
        [f_new,g_new,~] = feval(f_name,x_new);
        
        rho_k = (f - f_new)/(-gTd);
        
        if rho_k >= eta_s
            x = x_new;
            f = f_new;
            g = g_new;
            if rho_k >= eta_e
                delta_k = max([delta_k gamma_e*norm(d_k)]);
            end
        else
            delta_k = gamma_c*norm(d_k);
        end
        
        k = k + 1;

    end
    
    disp(toc(start_time))
    
    % display optimization parameters
    space = sprintf('\n');
    fprintf('\ntotal number of iterations:\n')
    disp(space)
    disp(k)
    fprintf('\nminimizing value of x:\n')
    disp(space)
    disp(x_new)
    fprintf('\nfunction minimum:\n')
    disp(space)
    disp(f_new)
    fprintf('\ngradient at minimizing value of x:\n')
    disp(space)
    disp(g_new)

end

