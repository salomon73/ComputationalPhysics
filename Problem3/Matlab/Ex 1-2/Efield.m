%% Electric field %%
function [Ex,Ey]=Efield(N,V)
    Ex=zeros(N,N); Ey=Ex; dx=N/(N+1); dx2=-1/(2*dx);
    for i=1:N
        for j=1:N
            if j~=1 && j~=N
                Ex(i,j)=dx2*(V(i,j+1)-V(i,j-1));
            elseif j==1 && i>=N/3+1 && i<=2*N/3
                    Ex(i,j)=2*dx2*(V(i,j+1)-V(i,j)); %backward
            elseif j==N && i>=N/3+1 && i<=2*N/3
                    Ex(i,j)=2*dx2*(V(i,j)-V(i,j-1)); %forward
            end
            if i~=1 && i~=N
                Ey(i,j)=dx2*(V(i+1,j)-V(i-1,j));
            end
        end
    end
end