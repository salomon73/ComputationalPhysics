function [L,U,P] = lu_decomposition(A)
sz=size(A);
L=eye(sz(1));
P=eye(sz(1));
U=A;
    if sz(1)~=sz(2)
        error('This is not a square matrix.');
    else
         for i=1:sz(1)
             [~,r] = max(abs(U(i:sz(1),i)));
             r=r+i-1;
             if r>=i
                 U([r,i],:)=U([i,r],:);
                 P([r,i],:)=P([i,r],:);
             end
             if i>=2
                L([r,i],1:i-1)=L([i,r],1:i-1);
             end
            for j=i:sz(1)-1
                L(j+1,i) = U(j+1,i)/U(i,i);
                U(j+1,:) = U(j+1,:)-L(j+1,i)*U(i,:); 
            end
        end
    end
end