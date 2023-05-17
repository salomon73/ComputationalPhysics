%% [VEC, VAL] = EIG_IPOWER(A,T)
% =============================================
% 
% Solves the eigenvalue problem Ax = lambda x
% for the eigenvalue lambda and the eigenvector x
% with an estimate of the eigenvalue t being 
% provided, with inverse power method
% 
% INPUT
% -----
%   -A:  2D complex & Hermitian matrix 
%        for the eigenvalue problem
%   -t:  estimated value for the eigenvalue lambda
%
% OUTPUT
% ------
%   -vec: eigenvector for lambda 
%   -val: eigenvalue lambda
% ------------------------------------%
% Written by S.Guinchard (05/03/23)   %
% ------------------------------------%
function [vec,val] = eig_ipower(A,t)
    I     = eye(size(A));
    shift = A-t*I;
    vec  = rand(length(A),1); % Initialize container for result
    vec  = vec./norm(vec);    % Normalize
    val  = ctranspose(vec)*A*vec;
    temp = val + 1;  

    while abs(val-temp)>eps 
        vec = inv(shift)*vec;
        vec = vec/norm(vec);
        temp = val;
        val = ctranspose(vec)*A*vec;       
    end
end