% An initialization function for finding an initial vertex of linear
% inequality constraints.
%
%   1. Calls modified simplex method as a subroutine.
%
%   2. Solves the problem
%
%             Ax >= b
%
%      with inputs [A,b,x] for any initial x.
%
%   3. Returns 'infeasible' if the feasible set is null.

% Jeremy Orosco
% MATH 271C, Spring 2016
% Coimbra Energy Group
% Department of Mechanical
% and Aerospace Engineering
% University of California, San Diego

% Copyright 2016 Jeremy Orosco
% This work is licensed under the
% Creative Commons Attribution-NonCommercial-ShareAlike 4.0
% International License. To view a copy of this license,
% visit http://creativecommons.org/licenses/by-nc-sa/4.0/.


function x_feas = phase_one(A,b,x)

    % get system sizes
    n = size(A,2);
    m = size(A,1);
    
    % define (approximate zero) tolerance
    tolerance = 10e01*sqrt(eps);
    
    % preallocate
    theta = zeros(m,1);
    f = zeros(m,1);
    
    % define initial feasible point if not given
    % (tries to initialize to a vertex)
    if nargin < 3
    
        % get indep. working set from active set
        [~,~,p] = lu(A,'vector');
        mw = rank(A);
        working = p(1:mw);
        Aw = A(working,:);
        bw = b(working);
        
        % compute init. point
        x = linsolve(Aw,bw);

    end
    
    % check feas. (i.e. vertex)
    if min(A*x-b) + tolerance >= 0
        x_feas = x;
    else
        % compute resid.
        r = A*x-b;
        
        % compute max constraint violation
        for i = 1:m
           theta(i) = max(-r(i),0);
           if theta(i) > tolerance
            f(i) = max(sign(-r(i)),0);
           end
        end
        theta0 = max(theta);
        
        % define augmented (phase-1) system
        A_ = [A          f;
              zeros(1,n) 1];
        b_ = [b;
              0];
        g_ = [zeros(n,1);
              1];
        x_ = [x;
              theta0];
          
        % run modified simplex method on augmented system
        [x_s,~] = modified_simplex(A_,b_,g_,x_);
        
        % if vertex found, return result
        if abs(x_s(end)) < tolerance
            x_feas = x_s(1:end-1);
        % else, no feasible point
        elseif abs(x_s(end)) >= tolerance
            x_feas = Inf;
        end
    end

end




