function [f,g,H] = rosenbrock(x)

a = 1;
b = 100;

f = (a-x1)^2+b*(x2-x1^2)^2;

sqrt((t.^(1-t)./gamma(2-t)-(2*sqrt(t))/sqrt(pi))^2)

g1 = 2*x1-2*a+4*b*x1^3-4*b*x1*x2;
g2 = -2*b*x1^2+2*b*x2;
g = [g1;g2];

h11 = 2+8*b*(x1^2)-4*b*(x2-x1^2);
h12 = -4*b*x1;
h21 = -4*b*x1;
h22 = 2*b;

H = [h11 h12;
     h21 h22];

end

