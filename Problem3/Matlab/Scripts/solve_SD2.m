function [x,Err]=solve_SD2(A,b,n)
% Arguments : A (square symmetric matrix), b (real vector), n (positive
% integer)
% Computes the solution x of Ax=b using the steepest descent method, and
% the relative error at each iteration
% Returns : x (real vector), Err (1D array)
    N = length(b);
    x = ones(N,1);
    r = b - A*x;
    % Exact solution
    X = A\b;
    % Relative error w.r.t. to X
    Err0 = sqrt((x-X)'*A*(x-X));
    Err = 1;
    k = 0;
    while k < n
        alpha = (r'*r)/(r'*A*r);
        x = x + alpha*r;
        r = b - A*x;
        % Relative error w.r.t. to X
        Err = [Err sqrt((x-X)'*A*(x-X))/Err0];
        k=k+1;
    end
end