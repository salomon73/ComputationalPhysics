%% Discretization of Poisson's eq %% 
function [A,rho]=discretization(N,L,V,C)
   A=zeros(N^2,N^2); rho=zeros(N,N);
   dx=L/(N+1); dx2=1/dx^2;
   %building matrix A
   for i=1:N^2
        if i==1
            A(i,i+1)=dx2;
            A(i,i+N)=dx2;
            A(i,i)=-4*dx2;
        elseif i==N
            A(i,i-1)=dx2;
            A(i,i+N)=dx2;
            A(i,i)=-4*dx2;
        elseif i==N*(N-1)+1
            A(i,i+1)=dx2;
            A(i,i-N)=dx2;
            A(i,i)=-4*dx2;
        elseif i==N^2
            A(i,i-1)=dx2;
            A(i,i-N)=dx2;
            A(i,i)=-4*dx2;
        elseif i<N
            A(i,i-1)=dx2;
            A(i,i+1)=dx2;
            A(i,i+N)=dx2;
            A(i,i)=-4*dx2;
        elseif i>N*(N-1)+1
            A(i,i-1)=dx2;
            A(i,i+1)=dx2;
            A(i,i-N)=dx2;
            A(i,i)=-4*dx2;
        elseif 1==mod(i,N)
            A(i,i+1)=dx2;
            A(i,i+N)=dx2;
            A(i,i-N)=dx2;
            A(i,i)=-4*dx2;
        elseif 0==mod(i,N)
            A(i,i-1)=dx2;
            A(i,i-N)=dx2;
            A(i,i+N)=dx2;
            A(i,i)=-4*dx2;
        else
            A(i,i-1)=dx2;
            A(i,i+1)=dx2;
            A(i,i-N)=dx2;
            A(i,i+N)=dx2;
            A(i,i)=-4*dx2;
        end
   end
   %there is a charge at the center of the grid
   for i=1:N^2
       if i==(N^2+1)/2
           A(i,i)=dx2;
           A(i-1,i)=0;
           A(i+1,i)=0;
           A(i+N,i)=0;
           A(i-N,i)=0;
       else
           A((N^2+1)/2,i)=0;
       end
   end
   %puting B.C in rho
   for i=N/3+1:2*N/3
       rho(1,i)=-V*dx2;
       rho(end,i)=-V*dx2;
   end
   %considering the B.C due to the charge in the center
   rho((N^2+1)/2)=C*dx2;
   rho((N^2+1)/2+1)=-C*dx2;
   rho((N^2+1)/2-1)=-C*dx2;
   rho((N^2+1)/2+N)=-C*dx2;
   rho((N^2+1)/2-N)=-C*dx2;
   %continuity between the conductor and the insulator
   rho(1,N/3+1)=-0.5*dx2; 
   rho(end,N/3+1)=-0.5*dx2;
   rho(1,2*N/3)=-0.5*dx2;
   rho(end,2*N/3)=-0.5*dx2;
   rho=rho';
   rho=reshape(rho,N^2,1);
end