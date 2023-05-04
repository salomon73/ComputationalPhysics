function [L,U] = lu_decomposition_nopiv(A)
sz=size(A);
L=eye(sz(1));
U=A;
    if sz(1)~=sz(2)
        error('This is not a square matrix.');
    else
        for i=1:sz(1)
            for j=i:sz(1)-1
                L(j+1,i) = U(j+1,i)/U(i,i);
                U(j+1,:) = U(j+1,:)-L(j+1,i)*U(i,:); 
            end
        end
    end
end