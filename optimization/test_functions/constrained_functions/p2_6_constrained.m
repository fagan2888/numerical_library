function [f,c,g,J,H,G,Hc] = p2_6_constrained(x,y)

    if nargin < 2
        y = [0;0];
    end

    x0 = [1;3;5];
    A = [2 7 1;
         2 3 3];
    b = [1;4];

    f = x'*x;
    g = 2*(x-x0);
    H = 2*eye(numel(x));

    c = A*x - b;
    J = A;

    G = [g - J'*y;
         -c       ];

    Hc = zeros(3,3,2);
  
end