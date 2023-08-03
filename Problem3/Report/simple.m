function x = solve_SD(A,b)
% Computes the solution x of Ax=b using the steepest descent method
% INPUT  : A (square symmetric matrix), b (1D vector)
% OUTPUT : x (1D vector)
    N = length(b);
    x = ones(N,1);
    r = b - A*x;
    while norm(r) > 100*eps
        alpha = (r'*r)/(r'*A*r);
        x = x + alpha*r;
        r = b - A*x;
    end
end