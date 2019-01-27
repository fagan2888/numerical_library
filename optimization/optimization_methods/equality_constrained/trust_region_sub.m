% Trust-region subproblem solver.
%
%   1. Solves constrained optimization trust-region subproblem
%
%           minimize g'*d + (1/2)*d'*H*d
%
%         subject to norm(d) <= delta
%
%      where delta is the trust-region radius.
%
%   2. Includes conditions for addressing the hard case.

% Jeremy Orosco
% MATH 271B, Fall 2016
% Coimbra Energy Group
% Department of Mechanical
% and Aerospace Engineering
% University of California, San Diego

% Copyright 2016 Jeremy Orosco
% This work is licensed under the
% Creative Commons Attribution-NonCommercial-ShareAlike 4.0
% International License. To view a copy of this license,
% visit http://creativecommons.org/licenses/by-nc-sa/4.0/.


function [dp,dn] = trust_region_sub(g,B,delta)

    % initialized auxiliary values
    n = numel(g);                   % system size
    tolerance = 10e01*sqrt(eps);    % approximate zero radius
    epsilon = 10e-06;               % start dist. from most neg. e.value
    dn = [];                        % initialize second solution empty set
    
    % define definiteness of Hessian
    [Vn,Sn] = eigs(B,1,'SA');
    [~,S1] = eigs(B,1,'LA');
    positive_definite = Sn > tolerance;
    indefinite = Sn < -tolerance && S1 > tolerance;
    
    % if pos. def. compute trivial step
    if positive_definite
        dn = -linsolve(B,g);
    end
    
    % check for trivial solution
    if positive_definite && norm(dn) <= delta
        dp = dn;
        return
        
    % if not trivial, check degeneracy
    elseif indefinite || (positive_definite && norm(dn) > delta)
            
        % nondegenerate case
        if abs(Vn'*g) > tolerance

            % find sigma -> Newton's method
            sigma_new = -Sn+epsilon;
            sigma = sigma_new + 1;
            while abs(sigma_new-sigma) > tolerance
                sigma = sigma_new;
                d_sigma = -linsolve(B+sigma*eye(n),g);
                uTu = d_sigma'*((B+sigma*eye(n))^-1)*d_sigma;
                f1 = (norm(d_sigma)^2)/uTu;
                f2 = (norm(d_sigma)-delta)/delta;
                sigma_new = sigma + f1*f2;

            end
            sigma = sigma_new;

            % return nondegenerate step
            dp = -linsolve(B+sigma*eye(n),g);
                
        % degenerate case
        elseif abs(Vn'*g) < tolerance

            % compute degenerate step
            d_dagger = -pinv(B-Sn*eye(n))*g;

            % easy case (degenerate step to boundary)
            if norm(d_dagger) >= delta
                dp = d_dagger;
                return

            % hard case (null component, multiple solutions)
            elseif norm(d_dagger) < delta

                % get unit null vector
                Z = null(B-Sn*eye(n));
                z = Z(:,1);
                
                % calculate null step magnitude
                tau = sqrt(delta^2-norm(d_dagger)^2);

                % return multiple solutions
                dp = d_dagger + abs(tau)*z;
                dn = d_dagger - abs(tau)*z;

            end

        end
            
    end
        
end