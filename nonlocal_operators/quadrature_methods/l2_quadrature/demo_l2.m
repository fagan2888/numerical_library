% A test script for the numerical evaluation of the fractional or
% variable order operator using an L2 quadrature.

% Jeremy Orosco
% Coimbra Energy Group
% Department of Mechanical
% and Aerospace Engineering
% University of California, San Diego

% define a temporal grid
T = 2*pi;
h = 0.05;
N = round(T/h);
t = (0:N)'*h;

% define derivative order
p = 0.499 - (0.5*sin(t)).^2;

% define function
A = 1;
w = 5;
f = A*cos(w*t);

% preallocate
[dp_actual,d1f,dp_approx] = deal(zeros(N,1));

% evaluate actual (looped due to poorly behaved hypergeometric function)
for n = 1:N+1
%     dp_actual(n) = ((t(n)^(1-p(n)))/gamma(2-p(n)))*hypergeom(1,[1-p(n)/2,3/2-p(n)/2],-(t(n)^2)/4);
    dp_actual(n) = c_sinusoid(t(n),p(n),A,w,'cos');
end

% start time
start_time = tic;

% get numerical approximation
for n = 1:N
    
    % create shifted index iterate
    ni = n + 1;
    
    % get fractional derivative
    dp_approx(ni) = voo_l2(n,h,p(ni),f(1:ni));
    
end

% display elapsed time
fprintf('elapsed time: %.2f', toc(start_time))

% compute the error and print to console
global_error = sum(abs(dp_approx-dp_actual))/N;
fprintf('\nglobal error: %.4e\n\n', global_error)

% plot results
figure(1)
clf
hold on
box on
grid on
plot(t,dp_actual,'k--','linewidth',3.5)
plot(t,dp_approx,'linewidth',2)
axis([0 t(end) -1.1*abs(min(dp_actual)) 1.1*abs(max(dp_actual))])
legend({'actual','approximate'},'fontsize',14,'location','north')

