%% Creation of EWV %%
function EWV=ewv(M)
[l,c]=size(M);
M=pinv(M);
EWV=0; temp=0;
    for i=1:c
        for j=1:l
            temp=temp+M(i,j)^2*(0.1+0.001*(j-1))^2;
        end
        EWV=EWV+temp;
    end
end