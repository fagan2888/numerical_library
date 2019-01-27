% Constrained example function:
%
%   f(x) = x1*x2 - (1/2)*x1^2-x2
%
% subject to
%
%   c(x) = 0
%
% with
%
%   c(x) = r^2 - x1^2 - x2^2
%
% suprema (approximate) for r = sqrt(2):
%
%   minimizer (global)  :	f(-1,1) = -2.5
%   minimizer (local)   :   f(1.38,-0.31) = -1.07
%   maximizer (global)  :   f(-0.58,-1.29) = 1.87
%   maximizer (local)   :   f(1,1) = -0.5
%
% returns evaluations:
%
%   f: objective function
%   c: constraint vector
%   g: gradient of objective function
%   J: Jacobian matrix of constraint gradients
%   H: Hessian of Lagrangian, w.r.t. x
%   G: gradient of Lagrangian, w.r.t. (x,y)
%   Hc: Hessian of constraint vector (in general, a 3-D matrix)

% Jeremy Orosco
% MATH 271B, Winter 2016
% Coimbra Energy Group
% Department of Mechanical
% and Aerospace Engineering
% University of California, San Diego


function [f,c,g,J,H,G,Hc] = example1_constrained(x,y)

    if nargin < 2
        y = 0;
    end

    x1 = x(1);
    x2 = x(2);

    r = sqrt(2);

    f = x1*x2-(1/2)*x1^2-x2;
    c = r^2 - x1^2 - x2^2;

    g1 = x2-x1;
    g2 = x1-1;
    g = [g1;
         g2];

    J1 = -2*x1;
    J2 = -2*x2;
    J = [J1 J2];

    H11 = -1;
    H12 = 1;
    H21 = 1;
    H22 = 0;
    H = [H11 H12;
         H21 H22];

    G = [g-J'*y;
         -c    ];

    Hc11 = -2;
    Hc12 = 0;
    Hc21 = 0;
    Hc22 = -2;
    Hc = zeros(2,2,1);
    Hc(:,:,1) = [Hc11 Hc12;
                 Hc21 Hc22];

end

