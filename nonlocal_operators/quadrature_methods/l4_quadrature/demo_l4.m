% A test script for the numerical evaluation of the fractional or
% variable order operator using an L4 quadrature
% as given in Cao et al. (2015).

% Jeremy Orosco
% Coimbra Energy Group
% Department of Mechanical
% and Aerospace Engineering
% University of California, San Diego

% choose variable or fractional order system
%
%   -> 'fractional' for constant fractional order
%   -> 'variable' for variable order
%
derivative_type = 'variable';

% choose which function to evaluate
%
%   -> 'f1' for t^4
%   -> 'f2' for sin(t)
%
which_function = 'f2';

% define a temporal grid
T = 2*pi;
h = 0.05;
N = round(T/h);
t = (0:N)'*h;

% get derivative vector
switch derivative_type
    case 'fractional'
        p = 0.65*ones(N+1,1);
    case 'variable'
        p = 0.499 - (0.5*sin(t)).^2;
    otherwise
        disp('must choose a valid derivative type...')
        return
end

% get function vector and analytical vector
switch which_function
    case 'f1'
        dp_approx = zeros(N,1);
        f = t.^4;
        dp_actual = (24./gamma(5-p)).*t.^(4-p);
    case 'f2'
        [dp_actual,dp_approx] = deal(zeros(N,1));
        f = sin(t);
        for n = 1:N+1
            dp_actual(n) = ((t(n)^(1-p(n)))/gamma(2-p(n)))*hypergeom(1,[1-p(n)/2,3/2-p(n)/2],-(t(n)^2)/4);
        end
    otherwise
        disp('must choose a valid function...')
        return
end

% start time
start_time = tic;

% get numerical approximation
for n = 1:N
    
    % create shifted index iterate
    ni = n+1;
    
    % get fractional derivative
    dp_approx(ni) = voo_l4(n,h,p(ni),f(1:ni));
    
end

% display elapsed time
fprintf('elapsed time: %.2f', toc(start_time))

% compute the error and print to console
global_error = sum(abs(dp_actual-dp_approx))/N;
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

