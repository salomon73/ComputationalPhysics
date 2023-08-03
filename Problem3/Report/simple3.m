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