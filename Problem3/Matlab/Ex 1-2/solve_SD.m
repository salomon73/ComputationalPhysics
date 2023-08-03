%% Steepest descend method %%
function [X,count,err]=solve_SD(A,b)
    X=ones(length(b),1); Xref=A\b;
    r=b-A*X; count=0;
    err=zeros(length(b),1);
    while norm(r)>1e-11 && count<=1e8
        alpha=r'*r/(r'*A*r);
        X=X+alpha*r;
        r=b-A*X;
        count=count+1;
        err(count)=norm(Xref-X)/norm(Xref);
    end
end