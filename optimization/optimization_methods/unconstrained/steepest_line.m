function [x_new,f_new,g_new] = steepest_line(f_name,x)

    k = 0;
    gamma_c = 1/2;
    eta_s = 1/4;
    tol = 10*sqrt(eps);
    iter_limit = 10e03;
    
    start_time = tic;

    [f,g,~] = feval(f_name,x);

    while norm(g) > tol && k < iter_limit
        
        j = 0;
        
        p = -g;
        gTp = g'*p;
        
        alpha_k = 1;
        x_new = x + alpha_k*p;
        
        [f_new,g_new,~] = feval(f_name,x_new);

        while f_new > f + alpha_k*eta_s*gTp

            alpha_k = gamma_c*alpha_k;
            x_new = x + alpha_k*p;
            [f_new,g_new,~] = feval(f_name,x_new);
            j = j + 1;

        end
        
        x = x_new;
        f = f_new;
        g = g_new;
        
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

