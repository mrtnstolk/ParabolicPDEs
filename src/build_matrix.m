function [A] = build_matrix(D, h, k, n)
%BUILD_MATRIX Construct matrix "A" to solve Crank-Nicolson model.
% 
%   D := Diffusion constant.
%   h := Step length for x.
%   k := Step length for t.
%   n := Size of the matrix.
%
%   A := Resulting tridiagonal matrix.

% Values for diagonal entries.
r = (D * k) / h^2;
R = r * ones(n, 1);

% Fill sparse matrix "A" with diagonal entries.
A = spdiags([R -2*R R], -1:1, n, n);

% Neumann conditions at endpoints "a" and "b".
A(1, 2) = 2*r; 
A(end, end-1) = 2*r;

end