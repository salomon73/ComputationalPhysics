%% [VEC, VAL] = EIG_RQ(A,T)
% =============================================
% 
% Solves the eigenvalue problem Ax = lambda x
% for the eigenvalue lambda and the eigenvector x
% with an estimate of the eigenvalue t being 
% provided, with Reyleigh Quotient
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
function [vec,val] = eig_rq(A,t)
   
    I     = eye(size(A));    % Identity
    chara = A - t * I;       % Characteristic matrix (cf charact. poly.)
    vec = rand(length(A),1); % Generate random vector
    vec = vec./norm(vec);    % Normalize it
    val = ctranspose(vec)*A*vec;
    temp = t;       

    while abs(val-temp)>eps  % machine precision epsilon
        vec = (eye(size(chara))/(chara))*vec; % Invert the matrix
        vec = vec/norm(vec);
        temp = val;
        val = ctranspose(vec)*A*vec;       
    end

end