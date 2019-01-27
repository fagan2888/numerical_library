% Compute an FOO or VOO using a fourth order quadrature.
%
% Given a derivative order, 'p', at time 'tn' and a function, 'f',
% returns the corresponding fractional or variable order derivative for
% that moment in time. Uses the method given in Cao et al. (2015).
%
%   Inputs:
%
%       n : outer (VODE) loop variable corresponding to time tn
%       h : uniform time step
%       p : derivative order
%       f : historical function values since departure from equilibrium
%
%   Output:
%
%      dp : fractional or variable order derivative
%
%   Notes:
%
%       (1) For an outer loop at the nth iteration, must have already
%           solved and stored f(n+1) before evaluating.
%
%   Example usage:
%
%       dp = voo_l4(n,h,p,f(1:n+1));

% Jeremy Orosco
% Coimbra Energy Group
% Department of Mechanical
% and Aerospace Engineering
% University of California, San Diego

function dp = voo_l4(n,h,p,f)

    % fractional parameters
    p1 = (1-p);
    p2 = (2-p);
    p3 = (3-p);

    % initialize and compute kernel sum
    for i = 1:n

        % historical spacing parameters
        s = n-i;
        sp = (s+1);

        % computes first interval
        if i == 1
            
            % auxiliary parameters
            s1 = p3*p2*(sp^p1);
            s2 = p3*p2*(s^p1);

            % rectangle update
            a = s1-s2;

            % initialize kernel sum
            kernel_sum = a*(f(i+1)-f(i));

        % computes second interval
        elseif i == 2
            
            % auxiliary parameters
            s1 = p3*p2*(sp^p1);
            s2 = p3*p2*(s^p1);
            s3 = p3*(sp^p2);
            s4 = p3*(s^p2);

            % trapezoidal update
            a = s1-s2;
            b = s3-s4-(s1+s2)/2;

            % update kernel sum
            kernel_sum = kernel_sum + (a+b)*f(i+1)-(a+2*b)*f(i)+b*f(i-1);

        % computes all remaining intervals
        elseif i > 2
            
            % auxiliary parameters
            s1 = p3*p2*(sp^p1);
            s2 = p3*p2*(s^p1);
            s3 = p3*(sp^p2);
            s4 = p3*(s^p2);
            s5 = sp^p3;
            s6 = s^p3;

            % auxiliary vector
            sv = [s1;s2;s3;s4;s5;s6];

            % auxiliary matrix
            W = [(1/3) -(11/6)  1 -2  1 -1
                 (1/2)    3    -2  5 -3  3
                 -1     -(3/2)  1 -4  3 -3
                 (1/6)   (1/3)  0  1 -1  1];

            % compute update weights
            w = W*sv;

            % update kernel sum
            kernel_sum = kernel_sum + w'*flipud(f(i-2:i+1));
            
        % error otherwise
        else
            disp('problem computing L4 quadrature...')
            return

        end

    end
    
    % compute operator
    dp = ((h^(-p))/gamma(4-p))*kernel_sum;

end