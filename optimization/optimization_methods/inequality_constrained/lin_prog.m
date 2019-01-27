% An arbitrary initial (not necessarily feasible) linear program solver.
%
%   1. Calls `phase_one' and `modified_simplex' as subroutines.
%
%   2. Solves the inequality constrained LP
%
%               minimize g'x
%             subject to Ax >= b
%
%      with inputs [A,b,g(,x)], where (optional) x is any initial point.
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


function [x_star, f_star] = lin_prog(A,b,g,x)

    % start timer
    start_time = tic;
    
    % get initial feasible point
    if nargin < 4
        x_feas = phase_one(A,b,zeros(numel(g),1));
    else
        x_feas = phase_one(A,b,x);
    end
    
    % return error if infeasible
    if x_feas == Inf
        error('constraints are infeasible')
    end
    
    % solve linear program
    [x_star,~,f_star] = modified_simplex(A,b,g,x_feas);
    
    % output results
    if f_star > -Inf
    
        % display elapsed time
        fprintf('\ntime elapsed: %0.4f\n',toc(start_time))

        % display initial vertex
        fprintf('\ninitial vertex:\n\n')
        disp(x_feas)

        % display solution
        fprintf('optimal vertex:\n\n')
        disp(x_star)

        % display optimum
        fprintf('optimum: %0.4f\n\n',g'*x_star)
        
    else
        
        % unbounded solution
        disp('solution unbounded')
        
    end

end

