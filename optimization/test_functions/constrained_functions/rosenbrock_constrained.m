% Constrained example function:
%
%   f(x) = (a-x1)^2+b*(x2-x1^2)^2
%
% subject to
%
%   c(x) = 0
%
% with
%
%   c(x) = r^2 - x1^2 - x2^2
%
% suprema (approximate):
%
%   if b > 0, minimizer (global)    :   f(a,a^2) = 0
%   if b < 0, maximizer (global)    :   f(a,a^2) = 0
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


function [f,c,g,J,H,G,Hc] = rosenbrock_constrained(x,y)

    if nargin < 2
        y = 0;
    end

    x1 = x(1);
    x2 = x(2);

    a = 1;
    b = 100;
    r = sqrt(2);

    f = (a-x1)^2+b*(x2-x1^2)^2;
    
    c = r^2 - x1^2 - x2^2;

    g1 = -2*(a-x1)-4*b*x1*(x2-x1^2);
    g2 = 2*b*(x2-x1^2);
    g = [g1;
         g2];

    J1 = -2*x1;
    J2 = -2*x2;
    J = [J1 J2];

    H11 = 2 + 8*b*x1^2 - 4*b*(x2-x1^2);
    H12 = -4*b*x1;
    H21 = -4*b*x1;
    H22 = 2*b;
    H = [H11 H12;
         H21 H22];

    G = [g-J'*y;
         -c     ];

    Hc11 = -2;
    Hc12 = 0;
    Hc21 = 0;
    Hc22 = -2;
    Hc = zeros(2,2,1);
    Hc(:,:,1) = [Hc11 Hc12;
                 Hc21 Hc22];

end

