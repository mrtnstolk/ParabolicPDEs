function [ x ] = tridiag_solver( A, b )
%TRIDIAG_SOLVER Solve tridiagonal system "Ax=b".
%   Use Thomas Algorithm to solve system in O(n) operations.
%   
%   A := Square (NxN) tridiagonal matrix.
%   b := RHS column vector (Nx1).
%
%   x := Solution column vector (Nx1).

n = size(A,1);
x = zeros(n,1);

A(1,2) = A(1,2) / A(1,1);
b(1) = b(1) / A(1,1);

for i = 2:(n-1)
    tmp = A(i,i) - A(i-1,i) * A(i,i-1);
    A(i,i+1) = A(i,i+1) / tmp;
    b(i) = (b(i) - b(i-1) * A(i,i-1)) / tmp;
end

x(n) = (b(n) - b(n-1) * A(n,n-1)) / (A(n,n) - A(n-1,n) * A(n,n-1));

for i = (n-1):(-1):1
    x(i) = b(i) - A(i,i+1) * x(i+1);
end

end

