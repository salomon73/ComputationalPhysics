%% Reconstruction of the density matrix %%
function D=density_matrix(sigma,rho)
D0=0.5*eye(2,2); D=0;
    for i=1:3
        D=D+rho(i)*sigma{1,i};
    end
D=D0+D;
end