%% Conjugate gradient method %%
function [X,count,err]=solve_CG(A,b)
    X=rand(length(b),1); Xref=A\b;
    d=b-A*X; r=d; count=0;
    err=zeros(length(b),1);
    while norm(r)>eps
        alpha=r'*r/(d'*A*d);
        X=X+alpha*d;
        temp_r=r;
        r=r-alpha*A*d;
        beta=r'*r/(temp_r'*temp_r);
        d=r+beta*d;
        count=count+1;
        err(count)=norm(Xref-X)/norm(Xref);
        %norm(r);
    end
end