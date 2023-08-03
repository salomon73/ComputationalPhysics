%% Non Linear Conjugate Gradient Method %%
function x = non_linCG(x)
    r=-grad(x);
    d=r; 
    while norm(r)>1e-5
        %it1=it1+1;
        alpha=-grad(x)'*d/(d'*myhess(x,d)*d);
        while true
            %it2=it2+1;
            ratio=d'*myhess(x,d)*d;
            if ratio<=0
                x=x+0.01*d;
            else
                alpha=-grad(x)'*d/(d'*myhess(x,d)*d);
                x=x+alpha*d;
            end
            if norm(alpha*d)<=5e-12
                break
            end
        end
        temp=r;
        r=-grad(x);
        beta=r'*r/(temp'*temp);
        d=r+beta*d;
        norm(r)
    end
end 