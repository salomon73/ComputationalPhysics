%% Projection matrix %%
function P=projection(N)
    P=zeros(3*N*(N-1)-1,3*N*(N-1));
    for i=1:3*N*(N-1)
        if i<=N*(N-1)-1
            P(i,[1 i+1])=[-1 1];
        elseif i>=N*(N-1) && i<=3*N*(N-1)-1
            P(i,i+1)=1; 
        end
    end
end