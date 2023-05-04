function [val, diago] = eig_j(A)
%% [VAL, DIAGO] = EIG_J(A)
% =============================================
% 
% Find the eigenvalues of matrix A and 
% gives both the ordered eigenvalues 
% and the diagonalized form of A, using 
% the Jacobi algorithm
% 
% INPUT
% -----
%   -A:  Matrix for the eigenvalue problem
%
% OUTPUT
% ------
%   -val: ordered vector of eigenvalues
%   -diago: diagonalized form of A
% ------------------------------------%
% Written by S.Guinchard (05/03/23)   %
% ------------------------------------%
    N = length(A(:,1));

    while(off(A)> 1e-6)
        test = 0; 
        for ii=1:N-1
            for jj = ii+1:N
                if (abs(A(ii,jj))>test)
                    p = ii;
                    q = jj;
                    test = abs(A(ii,jj));
                end
            end
        end

    
    if (A(p,q) ~=0)
        tau = (A(q,q)-A(p,p))/(2*A(p,q));
        if tau>=0
            t = sqrt(1+tau^2)-tau;
        else 
            t = -(sqrt(1+tau^2)-tau);
        end
        c = 1/sqrt(1+t^2);
        s = t*c;

    else 
        c = 1;
        s = 0;
    end

    Jac = eye(N);
    Jac(p,p) =  c;
    Jac(q,q) =  c;
    Jac(p,q) =  s;
    Jac(q,p) = -s;
    A = ctranspose(Jac)*A*Jac;
    end
    val = diag(A);
    diago = A;
end

function norm =  off(A)
% Compute the norm of the off-diagonal 
% terms of matrix A
%______________________________________
    N    = length(A(1,:));
    norm = 0;
    for ii = 1:N
        for jj = 1:N
            if ii==jj

            else 
                norm = norm + A(ii,jj)^2;
            end
        end
    end
    norm = sqrt(norm);
end