% A test script for the numerical evaluation of the fractional or
% variable order operator using the generalized trapezoidal quadrature
% as given in Soon et al. (2005).

% Jeremy Orosco
% Coimbra Energy Group
% Department of Mechanical
% and Aerospace Engineering
% University of California, San Diego

% define a temporal grid
T = 2*pi;
h = 0.01;
N = round(T/h);
t = (0:N)'*h;

% define derivative order
p = 0.499 - (0.5*sin(t)).^2;

% define function
f = sin(t);

% preallocate
[dp_actual,d1f,dp_approx] = deal(zeros(N,1));

% evaluate actual (looped due to poorly behaved hypergeometric function)
for n = 1:N+1
    dp_actual(n) = ((t(n)^(1-p(n)))/gamma(2-p(n)))*hypergeom(1,[1-p(n)/2,3/2-p(n)/2],-(t(n)^2)/4);
end

% get numerical approximation
for n = 0:N
    
    % create shifted index iterate
    ni = n + 1;
    
    % simulate first derivative calculation
    d1f(ni) = cos(t(ni));
    
    % get fractional derivative
    dp_approx(ni) = voo_trap(n,h,p(ni),d1f);
    
end

% compute the error and print to console
global_error = sum(abs(dp_approx-dp_actual))/N;
fprintf('\nglobal error: %.10f\n\n', global_error)

% plot results
figure(1)
clf
hold on
box on
grid on
plot(t,dp_actual,'k--','linewidth',2.5)
plot(t,dp_approx,'linewidth',2)
axis([0 t(end) -1.1*abs(min(dp_actual)) 1.1*abs(max(dp_actual))])
legend({'actual','approximate'})

