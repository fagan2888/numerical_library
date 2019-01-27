clear all
clc

% bounded solution exists
A = [  0    -4     7     4    -2
      -4    -1     4    -2     2
       0     1     0     4     0
       0    -1     0    -4     0
       0     1    -1     2    -9
      -6     1     0    -1     3
       2     0    -7    -1     0
       1     0     0     0     0
       0     1     0     0     0
       0     0     1     0     0
       0     0     0     1     0
       0     0     0     0     1 ];

b = [  5
       2
       0
      -1
     -10
      -3
      -5
       0
       0
       0
       0
       0 ];

% unbounded solution
% A = [  0    -4     7     4    -2
%       -4    -1     4    -2     2
%       -8    -2     8    -4     4];
% 
% b = [  5
%        2
%        4];

% infeasible constraints
% A = [  0    -4     7     4    -2
%       -4    -1     4    -2     2
%        1     0     0     0     0
%       -1     0     0     0     0];
% 
% b = [  5
%        2
%        1
%        1];

g = [ 10
       1
      -2
     -10
       5 ];

start_time = tic;
linprog(g,-A,-b)
disp(toc(start_time))

start_time = tic;
[x,f] = lin_solve(A,b,g);
disp(toc(start_time))

% matlab's solver returns infeasible primal, unbounded dual
% linprog(g,A,b)



