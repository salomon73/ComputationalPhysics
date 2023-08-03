%% Main for ex 10-11 %%
%% Graphical settings
fs=20;
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontsize', fs);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. Over-defined system of linear eq
clear,clc;
A=[3 2; 4 5; 1 1]; b=[1;4;1];
Penrose_inv=pinv(A); 
x=Penrose_inv*b;
%graphical verification 
xaxis=-0.3+x(1):1e-3:x(1)+0.3;
yaxis=-0.3+x(2):1e-3:x(2)+0.3;
r=zeros(length(xaxis),length(yaxis));
for i=1:length(xaxis)
    for j=1:length(yaxis)
        r(i,j)=norm(A*[xaxis(i),yaxis(j)]'-b);
    end
end
figure
c=contourf(xaxis,yaxis,r,50);
hold on
plot(x(1),x(2),'rX','MarkerSize',8,'LineWidth',1)
colorbar
legend('Residuals','$(-0.4118,1.1373)$')
title('$||A\vec x-\vec b||$')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5. Quantum state tomography
clear,clc;
%theoretical part
rho1=0.5*[1,1;1,1]; rho2=0.5*[3,-1;-1,1];
p1=load('ps1.csv');
p2=load('ps2.csv'); 
[V,Diag]=eig(rho1);
%p=[p1;p2];
N1=length(p1); N2=length(p2);
M0_1=0.5*ones(N1,1); M0_2=0.5*ones(N2,1); 
sigma_x=[0,1;1,0]; sigma_y=[0,-1i;1i,0]; sigma_z=[1,0;0,-1];
sigma={sigma_x,sigma_y,sigma_z};
%basis vectors
sz_up=[1;0]; sz_down=[0;1];
sx_up=1/sqrt(2)*(sz_up+sz_down);
sx_down=1/sqrt(2)*(sz_up-sz_down);
sy_up=1/sqrt(2)*(sz_up+1i*sz_down);
sy_down=1/sqrt(2)*(sz_up-1i*sz_down);
basis={sx_up,sx_down,sy_up,sy_down,sz_up,sz_down};
M=measurement_matrix(length(basis),sigma,E_operator(basis));
rank(M)
%now we solve the system M*rho=p-M0 with Penrose inverse
Pinv=pinv(M); rho_1=Pinv*(p1-M0_1); rho_2=Pinv*(p2-M0_2);
D1=density_matrix(sigma,rho_1); 
D2=density_matrix(sigma,rho_2);
[V1,Diag1]=eig(D1);
%verify we found the right state
verif=D1-V1(:,2)*V1(:,2)';
%taking only the 4 first projection operators
sigma_4={sigma_x,sigma_y};
basis_4={sx_up,sx_down,sy_up,sy_down};
M_4=measurement_matrix(length(basis_4),sigma_4,E_operator(basis_4));
rank(M_4)
% Pinv_4=pinv(M_4); rho_1_4=Pinv_4*(p1-M0_1); rho_2_4=Pinv_4*(p2-M0_2);
% D1_4=density_matrix(sigma_4,rho_1_4); 
% D2_4=density_matrix(sigma_4,rho_2_4);
[~,singular_values,~]=svd(M_4);
%% 5.3 QST with experimental constraints
clc; clear
%find in which case the density matrix can be reconstructed
phi=0; N=100; t=linspace(0,(N-1)*2*pi/N,N); 
%theta=linspace(0,pi,N);
theta=[0,pi/2,pi/3];
psi=cell(1,length(t));
sing_val=cell(1,length(theta));
EWV=zeros(length(theta),1);
for k=1:length(theta)
    for i=1:length(t)
        psi{i}=[cos(theta(k)/2); exp(-1i*N*t(i)/(N-1))*sin(theta(k)/2)];
    end
    E=E_operator(psi);
    sigma_x=[0,1;1,0]; sigma_y=[0,-1i;1i,0]; sigma_z=[1,0;0,-1];
    sigma={sigma_x,sigma_y,sigma_z};
    M=measurement_matrix(N,sigma,E);
    [~,val,~]=svd(M);
    sing_val{k}=diag(val);
    EWV(k)=ewv(M);
end
theta_mins=islocalmin(EWV);
%%
for i=1:length(theta_mins)
    if theta_mins(i)==1
        theta_min=theta(i);
    end
end
%%
figure
semilogy(theta,EWV,'.')
xlabel('$\theta$ [rad]')
ylabel('EWV')
%%
p3=load('ps3.csv');
M0=0.5*ones(length(p3),1);
Pinverse=pinv(M); rho_coeff=Pinverse*(p3-M0);
rho=density_matrix(sigma,rho_coeff);
%%
M1=[0.5 0.5; 0.5 0.5];
M2=0.5*[3 -1; -1 1];
[V,D]=eig(M1);


