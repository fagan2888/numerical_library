% Jeremy Orosco
% MATH 271B, Winter 2016
% Coimbra Energy Group
% Department of Mechanical
% and Aerospace Engineering
% University of California, San Diego


function [f,c,g,J,H,G,Hc] = p2_1b_constrained(x,y)

if nargin < 2
    y = 0;
end

x1 = x(1);
x2 = x(2);

f = exp(x1)*(4*x1^2+2*x2^2+4*x1*x2+2*x2+1);
c = 4*x1-x2-6;

g1 = exp(x1)*(1+4*x1^2+6*x2+2*x2^2+4*x1*(2+x2));
g2 = exp(x1)*(2+4*x1+4*x2);
g = [g1;
     g2];
 
J1 = 4;
J2 = -1;
J = [J1 J2];

H11 = exp(x1)*(9+4*x1^2+10*x2+2*x2^2+4*x1*(4+x2));
H12 = 2*exp(x1)*(3+2*x1+2*x2);
H21 = 2*exp(x1)*(3+2*x1+2*x2);
H22 = 4*exp(x1);
H = [H11 H12;
     H21 H22];

Hc = zeros(2,2,1);
Hc(:,:,1) = [0 0;
             0 0];
  
G1 = g-J'*y;
G2 = -c;
G = [G1;
     G2];

end

