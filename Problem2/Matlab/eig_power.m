%% [VEC, VAL] = EIG_POWER(A)
% =============================================
% 
% Solves the eigenvalue problem Ax = lambda x
% for the eigenvalue lambda and the eigenvector x
% with an estimate of the eigenvalue t being 
% provided, with the power method
% 
% INPUT
% -----
%   -A:  2D complex & Hermitian matrix 
%        for the eigenvalue problem
%
% OUTPUT
% ------
%   -vec: eigenvector for lambda 
%   -val: eigenvalue lambda
% ------------------------------------%
% Written by S.Guinchard (05/03/23)   %
% ------------------------------------%
function [vec,val] = eig_power(A)

    vec = rand(length(A),1); % Initialize container for result
    vec = vec./norm(vec);    % Normalize it
    val = ctranspose(vec)*A*vec;

    temp = val + 1;
    while abs(val-temp)>eps % machine precision epsilon
        vec = A*vec;
        vec = vec/norm(vec);
        temp = val;
        val = ctranspose(vec)*A*vec;       
    end
end