% Constrained example function.

% Jeremy Orosco
% MATH 271B, Winter 2016
% Coimbra Energy Group
% Department of Mechanical
% and Aerospace Engineering
% University of California, San Diego


function [f,c,g,J,H,G,Hc] = maratos(x,y)

    if nargin < 2
        y = 0;
    end

    x1 = x(1);
    x2 = x(2);

    f = -x1+2*(x1^2+x2^2-1);
    c = 1-(x1^2+x2^2);

    g1 = 4*x1-1;
    g2 = 4*x2;
    g = [g1;
         g2];

    J1 = -2*x1;
    J2 = -2*x2;
    J = [J1 J2];

    H11 = 4;
    H12 = 0;
    H21 = 0;
    H22 = 4;
    H = [H11 H12;
         H21 H22];

    G = [g-J'*y;
         -c    ];

    Hc11 = 2;
    Hc12 = 0;
    Hc21 = 0;
    Hc22 = 2;
    Hc = zeros(2,2,1);
    Hc(:,:,1) = [Hc11 Hc12;
                 Hc21 Hc22];

end

