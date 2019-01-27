% Compute an FOO or VOO using the generalized trapezoidal method.
%
% Given a derivative order, 'p', at time 'tn' and a vector of historical
% first derivatives, 'd1f', of a function, 'f', returns the corresponding
% fractional or variable order derivative for that moment in time. Uses
% the method given in Soon et al. (2005).
%
%   Inputs:
%
%       n : outer (VODE) loop variable corresponding to time tn
%       h : uniform time step
%       p : derivative order
%     d1f : historical first derivatives since departure from equilibrium
%
%   Output:
%
%      dp : fractional or variable order derivative
%
%   Dependencies:
%
%       (1) Calls the function 'weight.m' which must return the appropriate
%           weight at each iteration.
%
%   Notes:
%
%       (1) For an outer loop at the nth iteration, must have already
%           solved and stored d1f(n+1) before evaluating.
%
%   Example usage:
%
%       dp = voo_trap(n,h,p,d1x);

% Jeremy Orosco
% Coimbra Energy Group
% Department of Mechanical
% and Aerospace Engineering
% University of California, San Diego

function dp = voo_trap(n,h,p,d1f)

    % define shifted index limit
    ni = n + 1;
    
    % initialize and populate weight vector
    w = zeros(ni,1);
    for i = 0:n
        
        % define shifted index iterable
        ii = i + 1;
        
        % store weight
        w(ii) = weight(n-i,n,p);
        
    end
    
    % compute operator
    dp = ((h^(1-p))/gamma(3-p))*w'*d1f(1:ni);

end