%% Lennard-Jones potential %%
function E = energy(x)
    E=0;
    N=length(x)/3; %N particles in 3D
    x=reshape(x,[3,N])'; 
    %Lennard-Jones contribution with epsilon=sigma=1
    for i=1:N
        for j=1:i-1
            nrm=norm(x(i,:)-x(j,:));
            ULJ=4*((1/nrm)^12-(1/nrm)^6);
            E=E+ULJ;
        end
    end
end