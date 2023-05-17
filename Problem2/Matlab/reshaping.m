%% Reshaping x-vector into V,Ir,It %%
function [V,Ir,It]=reshaping(N,x)
Nv=N*(N-1); Nir=2*N*(N-1); Nit=3*N*(N-1);
   V=(reshape(x(1:Nv),[N-1,N]));
   Ir=(reshape(x(Nv+1:Nir),[N-1,N]));
   It=(reshape(x(Nir+1:Nit),[N-1,N]));
% figure
% contourf(V,10);
% c=colorbar;
% c.Label.String = '[V]';
% title('Voltage $V$');
% figure
% contourf(Ir,10); 
% c=colorbar;
% c.Label.String = '[A]';
% title('Current $I_r$')
% figure
% contourf(It,10); 
% c=colorbar;
% c.Label.String = '[A]';
% title('Current $I_t$')
end