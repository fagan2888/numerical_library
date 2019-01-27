% Order 3-p convergent Caputo discrete 
% fractional derivative implementation.

% temporal grid
h = 1;
T = 5;
N = round(T/h);
t = (0:N)'*h;

% derivative order
pf = 0.5*ones(N+1,1);

% function
% f = t.^4;
A = 1;
w = 2.5;
f = A*cos(w*t);

% analytical derivative
% dp_actual = (24./gamma(5-pf)).*t.^(4-pf);
dp_actual = zeros(N+1,1);
for i = 1:N+1
    dp_actual(i) = c_sinusoid(t(i),pf(i),A,w,'cos');
end

% evaluation looop
dp_approx = zeros(N+1,1);
for n = 1:N

    % fractional parameters
    p = pf(n+1);
    p1 = (1-p);
    p2 = (2-p);
    
    ft = f(1:n+1);

    % initialize and compute kernel sum
    for i = 1:n

        % historical spacing parameters
        s = n-i;
        sp = (s+1);

        % computes first interval
        if i == 1

            % auxiliary parameters
            s1 = p2*(sp^p1);
            s2 = p2*(s^p1);

            % rectangle update
            a = s1-s2;

            % initialize kernel sum
            kernel_sum = a*(ft(i+1)-ft(i));

        % computes remaining intervals
        elseif i > 1

            % auxiliary parameters
            s1 = p2*(sp^p1);
            s2 = p2*(s^p1);
            s3 = (sp^p2);
            s4 = (s^p2);

            % trapezoidal update
            a = s1-s2;
            b = s3-s4-(s1+s2)/2;

            % update kernel sum
            kernel_sum = kernel_sum + (a+b)*ft(i+1)-(a+2*b)*ft(i)+b*ft(i-1);

        % error otherwise
        else
            disp('problem computing L3 quadrature...')
            return

        end

    end

    % compute operator
    dp_approx(n+1) = ((h^(-p))/gamma(3-p))*kernel_sum;

end

% compute and print global error
global_error = sum(abs(dp_approx-dp_actual))/N;
fprintf('global error: %.4e\n',global_error)

% plot results
figure(1)
clf
hold on
box on
grid on
plot(t,dp_actual,'k--','linewidth',3.5)
plot(t,dp_approx,'linewidth',2)



