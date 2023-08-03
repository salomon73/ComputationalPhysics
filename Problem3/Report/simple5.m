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