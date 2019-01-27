%TQ       Factorization routine for a real m by n matrix A.
%         [T,Q] = tq(A) produces an upper-triangular matrix T
%         of the same dimension as A and an orthogonal matrix Q
%         so that AQ = T.

%         Written  by Philip Gill,     UCSD, Oct 1993.

function [T,Q] = tq(A)

% QR decomposition of the transpose.

[m,n]   = size(A);
mrev    = [ m:-1:1 ];
nrev    = [ n:-1:1 ];
[Q,T]   = qr(A(mrev,:)');

% Determine effective nullity

tol   = eps*norm(A,'fro');
if m > 1
   d = sum(abs(T'));
else
   d = abs(T');
end
nul  = find(d <= tol);
rank = n - length(nul);
if rank ~= m
   error('tq.m fails: matrix has deficient row rank')
end

Q = Q(:,nrev);
T = T(nrev,mrev)';


