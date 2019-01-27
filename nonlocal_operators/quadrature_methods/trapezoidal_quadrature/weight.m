% Compute the temporal weighting for the generalized trapezoidal method.
%
% Given a derivative order, 'p', at time 'tn' and the "iterate distance", 
% 's', from the current time, returns the historical temporal weighting
% as outlined in Soon et al. (2005). Here
%
%                               s = n - i,
%
% where i is the iterable history index and corresponds to tn = n*h, with
% 'tn' being the current loop time.
%
%   Inputs:
%
%       s : iterate distance as defined above
%       n : iterate limit as defined above
%       p : derivative order
%
%   Output:
%
%       w : historical first derivative weighting
%
%   Example usage:
%
%       w = weight(n-i,n,p);

% Jeremy Orosco
% Coimbra Energy Group
% Department of Mechanical
% and Aerospace Engineering
% University of California, San Diego

% normal method
function w = weight(s,n,p)

    % get appropriate weight, catch errors
    if s == n
        w = ((s-1)^(2-p)-(n^(1-p))*(s+p-2))*(n~=0);
    elseif (0 < s) && (s < n)
        w = (s-1)^(2-p)-2*s^(2-p)+(s+1)^(2-p);
    elseif s == 0
        w = 1;
    else
        disp('weight indexing error...')
    end
    

end