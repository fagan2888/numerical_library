% The classical simplex method.
%
%   1. Solves the inequality constrained linear program
%
%           minimize g'x
%       subject to Ax >= b
%
% 	2. Requires that the user defined x is a feasible point.
%
%   3. Returns `unbounded' if the solution is unbounded.

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


function [x_star,f_star] = simplex(A,b,g,x)

    % define system size
    m = numel(b);
    
    % define tolerance
    tolerance = 10e01*sqrt(eps);
    
    % initialization / preallocatation
    k = 0;
    alpha_hat = zeros(m,1);

    % define residual
    r = A*x-b;

    % find active constraints
    active = find(abs(r) <= tolerance);
    Aa = A(active,:);
    mw = rank(Aa);
        
    % define initial vertex working set
    if k == 0
        [~,~,p] = lu(Aa,'vector');
        working = active(p(1:mw));
    end
    Aw = A(working,:);
        
    % get minimum multiplier
    y = linsolve(Aw',g);
    [~,s] = min(y);
    
    % run simplex method
    while y(s) < -tolerance
        
        % compute search direction
        es = zeros(mw,1); es(s) = 1;
        pk = linsolve(Aw,es);
        
        % compute maximum feasible step
        for i = 1:m
           if A(i,:)*pk <= -tolerance
               alpha_hat(i) = r(i)/(-A(i,:)*pk);
           else
               alpha_hat(i) = Inf;
           end
        end
        [alpha,t] = min(alpha_hat);
        
        % check for unbounded solution
        if alpha == Inf
            disp('unbounded')
            x_star = inf*x;
            f_star = -inf;
            return
        end
        
        % replace departing constraint with blocking constraint
        Aw(s,:) = A(t,:);
        working(s) = t;
        
        % update iterate
        x = x + alpha*pk;
        
        % define residual
        r = A*x-b;
        
        % get minimum multiplier
        y = linsolve(Aw',g);
        [~,s] = min(y);

        % update counter
        k = k + 1;
        
    end
    x_star = x;
    f_star = g'*x_star;
        
end