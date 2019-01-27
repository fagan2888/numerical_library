% Inequality constrained example function.

% Jeremy Orosco
% MATH 271C, Spring 2016
% Coimbra Energy Group
% Department of Mechanical
% and Aerospace Engineering
% University of California, San Diego


function [f,c,g,J,H,G,Hc] = example6_ie(x,y)

    if nargin < 2
        y = [0;0];
    end

    x1 = x(1);
    x2 = x(2);

    f = 3*x1-x2+x1^2+(3/2)*x2^2;
    
    c1 = -3*(x1-2)^2+x2+3;
    c2 = 2*(x1-2)^2+x2+3;
    c = [c1;
         c2];

    g1 = 3+2*x1;
    g2 = -1+3*x2;
    g = [g1;
         g2];

    J = [12-6*x1 1;
         -8+4*x1 1];

    H = [2  0;
          0 3];

    G = [g-J'*y;
         -c    ];

    Hc = zeros(2,2,2);
    Hc(1,1,1) = -6;
    Hc(1,1,2) = 4;

end

