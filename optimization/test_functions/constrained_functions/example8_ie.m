% Inequality constrained example function.

% Jeremy Orosco
% MATH 271C, Spring 2016
% Coimbra Energy Group
% Department of Mechanical
% and Aerospace Engineering
% University of California, San Diego


function [f,c,g,J,H,G,Hc] = example8_ie(x,y)

    if nargin < 2
        y = [0;0;0];
    end

    x1 = x(1);
    x2 = x(2);

    f = (1/4)*x1^2+x2^2+3*x1*x2-3*x1-6*x2;
    
    c1 = (-5/2)*x1+x2+19/4;
    c2 = x1-1/2;
    c3 = x2-1/2;
    c = [c1;
         c2;
         c3];

    g1 = (x1/2)+3*x2-3;
    g2 = 3*x1+2*x2-6;
    g = [g1;
         g2];

    J = [(-5/2) 1;
         1      0
         0      1];

    H = [1/2  3;
          3   2];

    G = [g-J'*y;
         -c    ];

    Hc = zeros(2,2,3);

end

