% Inequality constrained example function.

% Jeremy Orosco
% MATH 271C, Spring 2016
% Coimbra Energy Group
% Department of Mechanical
% and Aerospace Engineering
% University of California, San Diego


function [f,c,g,J,H,G,Hc] = example5_ie(x,y)

    if nargin < 2
        y = [0;0;0];
    end

    x1 = x(1);
    x2 = x(2);

    f = x1*x2-(1/2)*x1^2-x2;
    
    c1 = -3*x1+x2+3;
    c2 = 2*x1+x2+3;
    c3 = -x2+2;
    c = [c1;c2;c3];

    g1 = x2-x1;
    g2 = x1-1;
    g = [g1;
         g2];

    J = [-3  1;
          2  1;
          0 -1];

    H11 = -1;
    H12 = 1;
    H21 = 1;
    H22 = 0;
    H = [H11 H12;
         H21 H22];

    G = [g-J'*y;
         -c    ];

    Hc = zeros(2,2,3);

end

