function [x, Res] = casestudies_gmres(A, b, tol, maxIter)
% GMRES_CUSTOM Implements the GMRES algorithm for both matrices and operator objects.
%   [x, Res] = GMRES_CUSTOM(A, b, tol, maxIter) solves the linear system A*x = b
%   where A can be a matrix or an object implementing 'mtimes'. 'tol' is the
%   tolerance for the residual norm, and 'maxIter' is the maximum number of
%   iterations. The function returns the solution x and the residual norms Res.

if nargin < 3
    tol = 1e-6;  % Default tolerance
end
if nargin < 4
    maxIter = min(100, length(b));  % Default maximum iterations
end

n = length(b);
x = zeros(n, 1);  % Initial guess is the zero vector
Res = [];

r0 = b - A*x;
beta = norm(r0);
V = zeros(n, maxIter+1);
H = zeros(maxIter+1, maxIter);
V(:,1) = r0 / beta;

for j = 1:maxIter
    w = A*V(:,j);
    for i = 1:j
        H(i,j) = w' * V(:,i);
        w = w - H(i,j) * V(:,i);
    end
    H(j+1,j) = norm(w);
    if H(j+1,j) < eps
        break;  % Break if w is effectively zero
    end
    V(:,j+1) = w / H(j+1,j);
    
    % Solve the least squares problem
    e1 = zeros(j+1, 1); e1(1) = 1;
    y = H(1:j+1, 1:j) \ (beta * e1);  % Minimize norm(beta*e1 - H*y)
    x = V(:, 1:j) * y;  % Update the solution
    
    % Compute and store the current residual norm
    currentRes = norm(b - A*x) / norm(b);
    Res = [Res, currentRes];
    
    if currentRes < tol
        break;  % Convergence achieved
    end
end

end
