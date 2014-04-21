function [] = runtests()
%RUNTESTS Execute set of tests for Crank-Nicolson model.
%   Each test should be written in a seperate subfunction
%   in this file.

for i = 1:5
    n = randi([100 500]); % Size of matrix.
    test_tridiag_build(n);
    test_tridiag_solver(n,1e-8);
end

end

function [success] = test_tridiag_build(n)
% Check if a matrix is tridiagonal.
success = true;

% Initialise data.
A = build_matrix(rand(), rand(), rand(), n);

% Check that all entries off the main, upper, and lower diags are 0.
for i = 1:n
    for j = 1:n
        if i < (j-1) || i > (j+1)
            if A(i,j) ~= 0;
                warning('"test_tridiag_build" failed.');
                success = false;
                return;
            end
        end
    end
end

end

function [success] = test_tridiag_solver(n,tol)
% Compare thomas algorithm to matlab's "\" function.
success = true;

% Initialise data.
A = 2*speye(n) - build_matrix(rand(), rand(), rand(), n); % 2I-A
b = rand(n,1); % TODO: write "build_bvector" function and use here.

% Time each method.
tic; x_real = A \ b; toc;
tic; x_test = tridiag_solver(A,b); toc;

if abs(x_real - x_test) > tol
    warning('"test_tridiag_solver" failed.');
    success = false;
    return;
end

end
