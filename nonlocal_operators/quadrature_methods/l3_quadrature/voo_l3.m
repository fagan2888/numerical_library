% Compute an FOO or VOO using the generalized L3 quadrature.
%
% Given a derivative order, 'p', at time 'tn' and a function, 'f',
% returns the corresponding fractional or variable order derivative for
% that moment in time. Uses an L3 quadrature.
%
%   Inputs:
%
%       n : outer (VODE) loop variable corresponding to time tn
%       h : uniform time step
%       p : derivative order
%       f : historical first derivatives since departure from equilibrium
%
%   Output:
%
%      dp : fractional or variable order derivative
%
%   Notes:
%
%       (1) For an outer loop at the nth iteration, must have already
%           solved and stored d1f(n+1) before evaluating.
%       (2) The historical function vector 'f' should have size n+1. In
%           other words, must have already solved for f(n+1) in outer
%           method.
%
%   Example usage:
%
%       dp = voo_l3(n,h,p,f(1:n+1));

% Jeremy Orosco
% Coimbra Energy Group
% Department of Mechanical
% and Aerospace Engineering
% University of California, San Diego

function dp = voo_l3(n,h,p,f)

    % fractional parameters
    p1 = (1-p);
    p2 = (2-p);

    % initialize and compute kernel sum
    for i = 1:n

        % historical spacing parameters
        s = n-i;
        sp = (s+1);

        % computes first interval
        if i == 1
            
            % auxiliary parameters
            s1 = p2*(sp^p1);
            s2 = p2*(s^p1);

            % rectangle update
            a = s1-s2;

            % initialize kernel sum
            kernel_sum = a*(f(i+1)-f(i));

        % computes remaining intervals
        elseif i > 1
            
            % auxiliary parameters
            s1 = p2*(sp^p1);
            s2 = p2*(s^p1);
            s3 = (sp^p2);
            s4 = (s^p2);

            % trapezoidal update
            a = s1-s2;
            b = s3-s4-(s1+s2)/2;

            % update kernel sum
            kernel_sum = kernel_sum + (a+b)*f(i+1)-(a+2*b)*f(i)+b*f(i-1);
            
        % error otherwise
        else
            disp('problem computing L3 quadrature...')
            return

        end

    end
    
    % compute operator
    dp = ((h^(-p))/gamma(3-p))*kernel_sum;

end