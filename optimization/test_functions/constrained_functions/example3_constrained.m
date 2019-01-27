% Constrained example function

% Jeremy Orosco
% MATH 271B, Winter 2016
% Coimbra Energy Group
% Department of Mechanical
% and Aerospace Engineering
% University of California, San Diego


function [f,c,g,J,H,G,Hc] = example3_constrained(x,y)

    if nargin < 2
        y = [0;0;0];
    end

    x1 = x(1);
    x2 = x(2);
    x3 = x(3);
    x4 = x(4);
    x5 = x(5);

    f = (x1-1)^2 + (x1-x2)^2 + (x2-x3)^2 + (x3-x4)^2 + (x4-x5)^2;
    
    c1 = x1+x2^2+x3^3-2-3*sqrt(2);
    c2 = x2-x3^2+x4+2-2*sqrt(2);
    c3 = x1*x5-2;
    c = [c1;c2;c3];

    g1 =  4*x1-2*(1+x2);
    g2 = -2*(x1-2*x2+x3);
    g3 = -2*(x2-2*x3+x4);
    g4 = -2*(x3-2*x4+x5);
    g5 = -2*x4+2*x5;
    g  = [g1;
          g2;
          g3;
          g4;
          g5];

    J1 = [ 1, 2*x2,  3*x3^2, 0,  0];
    J2 = [ 0,    1, -2*x3,   1,  0];
    J3 = [x5,    0,     0,   0, x1];
    J = [J1;
         J2;
         J3];

    H = [ 4, -2,  0,  0,  0;
         -2,  4, -2,  0,  0;
          0, -2,  4, -2,  0;
          0,  0, -2,  4, -2;
          0,  0,  0, -2,  2];

    G = [g-J'*y;
         -c    ];

    
    Hc(:,:,1) = zeros(5);Hc(2,2,1)=2;Hc(3,3,1)=6*x3;
    Hc(:,:,2) = zeros(5);Hc(3,3,2)=-2;
    Hc(:,:,3) = zeros(5);[Hc(1,5,3),Hc(5,1,3)]=deal(1);

end

