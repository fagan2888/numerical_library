

% Jeremy Orosco
% Coimbra Research Group
% Department of Mechanical
% and Aerospace Engineering
% University of California, San Diego

% Copyright 2017 Jeremy Orosco
% This work is licensed under the
% Creative Commons Attribution-NonCommercial-ShareAlike 4.0
% International License. To view a copy of this license,
% visit http://creativecommons.org/licenses/by-nc-sa/4.0/.


function [x_new,f_new,g_new] = pdpf_line(f_name,x,y)

    % define optimization parameters
    nf = 0;                 % init. func. eval. counter
    gamma_c = 1/2;          % line-search contraction factor
    eta_s = 1/4;            % sufficient decrease ratio
    gamma = 10;             % barrier decrease factor
    mu = 1;                 % barrier parameter
    opts.SYM = true;        % use the symmetric solver on the Newton system
    itns_to_print = 200;    % number of iterations to print

    % get initial values
    [f,g,H,c,J,Hc] = feval(f_name,x);nf=nf+1;
    
    % check strict feasibility of initial point
    if min(c) <= 0
        fprintf('\nneed a strictly feasible initial point...\n\n')
        [x_new,f_new,g_new] = deal(0);
        return
    end
    
    % get system sizes
    n = numel(x);
    m = numel(c);

    % initialize multiplier estimate
    if (nargin < 3)
        y = 10e-03*ones(m,1);
    end
    
    % define limits
    global_tolerance = 10e01*sqrt(eps);
    local_tolerance = 10e03*sqrt(eps);
    [inner_limit,outer_limit,line_limit,feas_limit] = deal(100);
    e = ones(m,1);
    
    % outer path-following loop
    % converges when we arrive at constraint
    % with a minimized merit gradient
    [k,k_outer] = deal(0);
    while (norm(c.*y-mu*e) > global_tolerance) && (k_outer < outer_limit)
            
        % show trajectory progress
        if k > 0
                fprintf('----------------------------------------------------------------------------\n')
                disp('               moving to next point on primal-dual trajectory               ')
                fprintf('----------------------------------------------------------------------------\n')
        end
    
        % construct merit gradient
        C = diag(c);
        Y = diag(y);
        ye = mu*C^-1*e;
        Db = Y^-1*C;
        Gb = [g-J'*(2*ye-y);
              Db*(y-ye)];
          
        % inner merit minimization loop
        % converges when merit gradient norm is minimized
        % at current point in primal-dual trajectory
        k_inner = 0;
        while (norm(Gb) >= local_tolerance) && (k_inner < inner_limit)
        
            % print header iteration
            if k == 0
                fprintf('\n')
                fprintf(' itn   nf    step    distance     optimal        mu     merit       cost\n')
                fprintf('--------------------------------------------------------------------------\n')
                str1 = sprintf (' %3g %4g       - ', k, nf);
                if norm(c.*y-mu*e) < local_tolerance && norm(Gb) < local_tolerance
                    str2 = sprintf ('   (%7.1e)   (%7.1e)    %7.1e    %7s    %8.1e', norm(c.*y-mu*e), norm(Gb), mu, '--', f);
                elseif norm(c.*y-mu*e) < local_tolerance && norm(Gb) >= local_tolerance
                    str2 = sprintf ('   (%7.1e)    %7.1e   %7.1e    %7s    %8.1e', norm(c.*y-mu*e), norm(Gb), mu, '--', f);
                elseif norm(c.*y-mu*e) >= local_tolerance && norm(Gb) < local_tolerance
                    str2 = sprintf ('    %7.1e   (%7.1e)     %7.1e    %7s    %8.1e', norm(c.*y-mu*e), norm(Gb), mu, '--', f);
                elseif norm(c.*y-mu*e) >= local_tolerance && norm(Gb) >= local_tolerance
                    str2 = sprintf ('    %7.1e     %7.1e   %7.1e    %7s    %8.1e', norm(c.*y-mu*e), norm(Gb), mu, '--', f);
                end
                disp([str1 str2])
            end
            
            % build merit Hessian approximation
            % reflect and elongate eigenvalues as necessary
            % (i.e. generate positive definit approximation)
            E = mod_ldl(H);
            B = H + E;
            Hb = B - threeD_mult(y,Hc);
            
            % construct and solve path-following equations
            % use symmetric solver
            Kb = [Hb  J';
                  J -Db];
            rb = [g-J'*y;
                  Db*(y-ye)];
            v = linsolve(Kb,-rb,opts);
            
            % extract update directions
            % for primal-dual path following
            p = v(1:n);
            q = -v(n+1:n+m);

            % construct initial merit
            b1 = mu*sum(log(c));
            b2 = mu*(sum(log((y.*c)/mu)) + sum(e-(y.*c)/mu));
            fb = f - b1 - b2;
            
            % initialize alpha and find the largest feasible step
            alpha = 1;
            x_new = x + alpha*p;
            y_new = y + alpha*q;
            [f_new,g_new,H_new,c_new,J_new,Hc_new] = feval(f_name,x_new);nf=nf+1;
            k_feas = 0;
            while ((min(c_new) < 0) || (min(y_new) < 0)) && (k_feas < feas_limit)
                
                % decrease step size a little slower than line search
                alpha = alpha/1.5;
                x_new = x + alpha*p;
                y_new = y + alpha*q;
                
                % compute new objective
                [f_new,g_new,H_new,c_new,J_new,Hc_new] = feval(f_name,x_new);nf=nf+1;
                
                % break on error
                if k_feas == feas_limit
                    x_new = x;
                    disp('no feasible step')
                    break
                end
                
                % update the counter
                k_feas = k_feas + 1;
                
            end
            
            % construct initial merit update
            b1_new = mu*sum(log(c_new));
            b2_new = mu*(sum(log((y_new.*c_new)/mu)) + sum(e-(y_new.*c_new)/mu));
            fb_new = f_new - b1_new - b2_new;

            % perform line search if necessary
            k_line = 0;
            while (fb_new > fb + eta_s*alpha*(Gb'*[p;q]))
                
                % update step
                alpha = gamma_c*alpha;
                x_new = x + alpha*p;
                y_new = y + alpha*q;
                
                % compute new objective
                [f_new,g_new,H_new,c_new,J_new,Hc_new] = feval(f_name,x_new);nf=nf+1;
                
                % compute new merit
                b1_new = mu*sum(log(c_new));
                b2_new = mu*(sum(log((y_new.*c_new)/mu)) + sum(e-(y_new.*c_new)/mu));
                fb_new = f_new - b1_new - b2_new;
                
                % break on error
                if k_line == line_limit
                    x_new = x;
                    disp('no step taken')
                    break
                end
                
                % update line search loop counter
                k_line = k_line + 1;

            end
            
            % update the space
            x = x_new;
            y = y_new;
            f = f_new;
            g = g_new;
            H = H_new;
            c = c_new;
            J = J_new;
            Hc = Hc_new;
            C = diag(c);
            Y = diag(y);
            ye = mu*C^-1*e;
            Db = Y^-1*C;
            
            % compute the new merit gradient
            Gb_new = [g-J'*(2*ye-y);
                      Db*(y-ye)];
            
            % update convergence parameter
            Gb = Gb_new;

            % update inner and total loop counters
            k_inner = k_inner + 1;
            k = k + 1;
        
            % print remaining iterations
            if k <= itns_to_print
                str1 = sprintf (' %3g %4g %7.2f ', k, nf, alpha);
                if norm(c.*y-mu*e) < local_tolerance && norm(Gb) < local_tolerance
                    str2 = sprintf ('   (%7.1e)   (%7.1e)  %7.1e   %8.1e    %8.1e', norm(c.*y-mu*e), norm(Gb), mu, fb, f);
                elseif norm(c.*y-mu*e) < local_tolerance && norm(Gb) >= local_tolerance
                    str2 = sprintf ('   (%7.1e)    %7.1e   %7.1e   %8.1e    %8.1e', norm(c.*y-mu*e), norm(Gb), mu, fb, f);
                elseif norm(c.*y-mu*e) >= local_tolerance && norm(Gb) < local_tolerance
                    str2 = sprintf ('    %7.1e   (%7.1e)   %7.1e   %8.1e    %8.1e', norm(c.*y-mu*e), norm(Gb), mu, fb, f);
                elseif norm(c.*y-mu*e) >= local_tolerance && norm(Gb) >= local_tolerance
                    str2 = sprintf ('    %7.1e     %7.1e   %7.1e   %8.1e    %8.1e', norm(c.*y-mu*e), norm(Gb), mu, fb, f);
                end
                disp([str1 str2])
            end
            
        end

        % update barrier parameter
        mu = mu/gamma;
        
        % update outer loop counter
        k_outer = k_outer + 1;
    
    end

end

