% A simplex method modified for use with any initial feasible point.
%
%   1. Solves the inequality constrained LP
%
%               minimize g'x
%             subject to Ax >= b
%
%      with inputs [A,b,g,x], where x is any feasible point.
%
%   2. For any initial feasible point, generates an initial vertex by
%      defining n-mw linearly independent temporary constraints, where
%      mw < n is the rank of the active set matrix. If n = mw, the initial
%      feasible point is a vertex and the classical simplex method ensues.
%
%   3. Returns `infeasible' if the feasible set is null.
%
%   4. Returns `unbounded' if the objective decreases in an unbounded
%      feasible direction.

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


function [x,y,f] = modified_simplex(A,b,g,x)

    % get system size
    [m,n] = size(A);

    % initialize / preallocate
    tolerance = 10e02*sqrt(eps);
    alpha_hat = zeros(m,1);
    not_optimal = 1;
    
    % compute objective
    f = g'*x;

    % compute residual
    r = A*x-b;

    % define active set matrix
    active = abs(r) < tolerance;
    Aa = A(active,:);
    mw = rank(Aa);
    
    % identify linearly dependent columns (=> basis dimensions)
    [~,~,pc] = qr(Aa,'vector');
    pt = pc(mw+1:n);
    
    % define augmented constraints (=> add indep. temp. bounds corresp.
    % to dependent basis dimensions)
    I = eye(n);
    A_temp = [A;I(pt,:)];
    b_temp = [b;x(pt)];
    
    % define temporary bound indices
    temp_indices = (m+1:m+(n-mw));
    
    % compute augmented residual
    rt = A_temp*x-b_temp;
    
    % define augmented active set matrix
    active = find(abs(rt) < tolerance);
    Aa = A_temp(active,:);
    
    % define augmented working set matrix
    [~,~,pr] = lu(Aa,'vector');
    working_set = active(pr(1:n));
    Aw = A_temp(working_set,:);
        
    % define current temp. bound set
    [~,temp_current,~] = intersect(working_set,temp_indices);
    
    % compute multipliers
    y = linsolve(Aw',g);
    
    % if any, check optimality temp. bounds
    if numel(temp_current) > 0
        [yks,rs] = max(abs((y(temp_current))));
        s = temp_current(rs);
        if abs(yks) < tolerance
            not_optimal = 0;
        end
        
    % otherwise, check optimality of regular working set
    else
        [yks,s] = min(y);
        if yks > -tolerance
            not_optimal = 0;
        end
    end
    
    % modified simplex iterations
    while not_optimal
        
        % compute search direction
        es = zeros(n,1);es(s) = -sign(y(s));
        pk = linsolve(Aw,es);

        % compute residual
        r = A*x-b;

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
            x = x+Inf*pk;
            f = -Inf;
            break
        end
        
        % replace departing constraint with blocking constraint
        Aw(s,:) = A(t,:);
        working_set(s) = t;
        
        % update current temp. bound set
        [~,temp_current,~] = intersect(working_set,temp_indices);

        % update iterate and objective
        x = x + alpha*pk;
        f = g'*x;

        % update multipliers
        y = linsolve(Aw',g);
        
        % if any, check optimality temp. bounds
        if numel(temp_current) > 0
            [yks,rs] = max(abs((y(temp_current))));
            s = temp_current(rs);
            if abs(yks) < tolerance
                not_optimal = 0;
            end
        
        % otherwise, check optimality of regular working set
        else
            [yks,s] = min(y);
            if yks > -tolerance
                not_optimal = 0;
            end
        end
    
    end
    
end




