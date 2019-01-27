% Constrained example function:
%
%   f(x) = exp(x1) + x1^2 + 2*x1*x2 + 4*x2^2
%
% subject to
%
%       c(x) = 0
%
%   with
%
%       c(x) = x1 + 2*x2 - 6

% Jeremy Orosco
% MATH 271B, Winter 2016
% Coimbra Energy Group
% Department of Mechanical
% and Aerospace Engineering
% University of California, San Diego


function [f,c,g,J,H,G,Hc] = example2_constrained(x,y)

    if nargin < 2
        y = 0;
    end

    x1 = x(1);
    x2 = x(2);

    f = exp(x1)+x1^2+2*x1*x2+4*x2^2;
    c = (x1+2*x2-6);

    g1 = exp(x1)+2*x1+2*x2;
    g2 = 2*x1+8*x2;
    g = [g1;
         g2];

    J1 = 1;
    J2 = 2;
    J = [J1 J2];

    H11 = 2+exp(x1);
    H12 = 2;
    H21 = 2;
    H22 = 8;
    H = [H11 H12;
         H21 H22];

    G = [g-J'*y;
         -c    ];

    Hc11 = 0;
    Hc12 = 0;
    Hc21 = 0;
    Hc22 = 0;
    Hc = zeros(2,2,1);
    Hc(:,:,1) = [Hc11 Hc12;
                 Hc21 Hc22];

end

