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