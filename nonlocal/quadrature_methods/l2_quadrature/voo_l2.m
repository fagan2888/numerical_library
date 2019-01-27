% Compute an FOO or VOO using an L2 quadrature method.
%
% Given a derivative order, 'p', at time 'tn' and a function, 'f', returns
% the corresponding fractional or variable order derivative for that
% moment in time.
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
%       dp = voo_l2(n,h,p,f(1:n+1));

% Jeremy Orosco
% Coimbra Energy Group
% Department of Mechanical
% and Aerospace Engineering
% University of California, San Diego

function dp = voo_l2(n,h,p,f)
    
    % initialize and populate weight vector
    kernel_sum = 0;
    for i = 1:n
        
        % historical spacing
        s = n-i;
        
        % temporal weight
        w = (s+1)^(1-p)-s^(1-p);
        
        % update kernel sum
        kernel_sum = kernel_sum + w*(f(i+1)-f(i));
        
    end
    
    % compute operator
    dp = ((h^(-p))/gamma(2-p))*kernel_sum;

end