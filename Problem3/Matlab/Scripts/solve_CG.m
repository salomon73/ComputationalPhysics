function x=solve_CG(A,b)
% Computes the solution x of Ax=b using the conjugate method
% INPUT  : A (square symmetric matrix), b (1D vector)
% OUTPUT : x (1D vector)
    N = length(b);
    x = ones(N,1);
    r = b - A*x;
    d = r;
    u = 0;
    while norm(r) > 100*eps
        alpha = (r'*r)/(d'*A*d);
        x = x + alpha*d;
        r_ = r;
        r = r - alpha*A*d;
        beta = (r'*r)/(r_'*r_);
        d = r + beta*d;
        u = u+1;
    end
    disp(u)
end