%% Verify Kirschoff's law %%
% Plot tje current showing 
% same value flowing in one node
% and out of the other

function verif_kirschoff(Ir, It)
temp=zeros(size(Ir));
[r,c]=size(Ir);
    for i=1:r
        for j=1:c
            temp(i,j)=Ir(i,j)+It(i,j);
        end
    end
    for i=1:r
        for j=2:c
            temp(i,j)=temp(i,j)-Ir(i,j-1);
        end
    end
    for i=2:r
        for j=1:c
            temp(i,j)=temp(i,j)-It(i-1,j);
        end
    end
%figure
    figure
    contourf(temp); 
    c=colorbar;
    c.Label.String = '[A]';
end