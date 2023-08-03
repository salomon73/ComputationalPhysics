function [x,Err]=solve_CG2(A,b,n)
% Arguments : A (square symmetric matrix), b (real vector), n (positive
% integer)
% Computes the solution x of Ax=b using the steepest descent method, and
% the relative error at each iteration
% Returns : x (real vector), Err (1D array)
    N = length(b);
    x = ones(N,1);
    r = b - A*x;
    d = r;
    % Exact solution
    X = A\b;
    % Relative error w.r.t. to X
    Err0 = sqrt((x-X)'*A*(x-X));
    Err = 1;
    k = 0;
    while k < n
        alpha = (r'*r)/(d'*A*d);
        x = x + alpha*d;
        r_ = r;
        r = r - alpha*A*d;
        beta = (r'*r)/(r_'*r_);
        d = r + beta*d;
        % Relative error w.r.t. to X
        Err = [Err sqrt((x-X)'*A*(x-X))/Err0];
        k=k+1;
    end
end