%% Creating the measurement matrix %%
function M=measurement_matrix(m,sigma,E)
%here m is the # of basis vector
    M=zeros(m,3);
    for i=1:m
        for j=1:length(sigma)
            M(i,j)=trace(sigma{1,j}*E{1,i});
        end
    end
end