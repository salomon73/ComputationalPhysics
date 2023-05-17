function [val, diago] = eig_cj(A)
%% [VAL, DIAGO] = EIG_CJ(A)
% =============================================
% 
% Find the eigenvalues of matrix A and 
% gives both the ordered eigenvalues 
% and the diagonalized form of A, using 
% the cyclic Jacobi algorithm
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
    N = length(A(1,:));
    while(off(A) > 1e-6)
        for ii = 1:N-1
            for jj = ii+1:N
                if (A(ii,jj) ~=0)
                   tau = (A(jj,jj)-A(ii,ii))/2/A(ii,jj);
                   if tau>=0
                       t = sqrt(1+tau^2)-tau;
                   else 
                       t = -(sqrt(1+tau^2)+tau);
                   end
                   c = 1/sqrt(1+t^2);
                   s = t*c;
                else 
                    c = 1;
                    s = 0;
                end
                Jac = eye(N); % Builds the eigenbasis change matrix
                Jac(ii,ii) =  c;
                Jac(jj,jj) =  c;
                Jac(ii,jj) =  s;
                Jac(jj,ii) = -s;
                A = ctranspose(Jac)*A*Jac;
            end
        end
    end
    val   = diag(A); % take diagonal elements = eigenvalues 
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