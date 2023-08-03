%% Main for ex 10-11 %%
%% Graphical settings
fs=20;
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontsize', fs);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Solving systems of lin. eq. 

[x1_SD, count1_SD, error1_SD]=solve_SD(A1,b1);
[x1_CG, count1_CG, error1_CG]=solve_CG(A1,b1);
[x2_SD,count2_SD, error2_SD]=solve_SD(A2,b2);
[x2_CG,count2_CG, error2_CG]=solve_CG(A2,b2);

N1_SD=1:1:count1_SD; N2_SD=1:1:count2_SD;
N1_CG=1:1:count1_CG; N2_CG=1:1:count2_CG;
%figures
figure
loglog(N1_SD,error1_SD,'r-');
xlabel('Number of iterations $N_{it}$')
ylabel('$||x_A-x_M||/||x_M||$')
grid on
hold on
loglog(N1_CG,error1_CG,'b-');
xlabel('Number of iterations $N_{it}$')
ylabel('$||x_A-x_M||/||x_M||$')
legend('Steepest descent','Conjugate gradient')
title('$N=5$')
figure
loglog(N2_SD,error2_SD,'r-');
xlabel('Number of iterations $N_{it}$')
ylabel('$||x_A-x_M||/||x_M||$')
grid on
hold on
loglog(N2_CG,error2_CG,'b-');
xlabel('Number of iterations $N_{it}$')
ylabel('$||x_A-x_M||/||x_M||$')
legend('Steepest descent','Conjugate gradient')
title('$N=50$')
%reaching machine precision with SD method
kappa1=cond(A1); kappa2=cond(A2);
N_max1=log(eps)/log((kappa1-1)/(kappa1+1));
N_max2=log(eps)/log((kappa2-1)/(kappa2+1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Poisson's equation
clear,clc;
N=45; L=N; V=1; C=1; eps0=1;
[A,rho]=discretization(N,L,V,C); %input : N,L,V,C,eps0 
%solving
[x_SD_Poisson,countSD_Poisson,errorSD_Poisson]=solve_SD(A,rho);
[x_CG_Poisson,countCG_Poisson,errorCG_Poisson]=solve_CG(A,rho);
x_ref=A\rho; 
diffSD=norm(x_ref-x_SD_Poisson)/norm(x_ref);
diffCG=norm(x_ref-x_CG_Poisson)/norm(x_ref);
%ploting
x_ref_resh=reshaping(x_ref,N,'Matlab solution');
x_SD_resh=reshaping(x_SD_Poisson,N,'SD solution');
x_CG_resh=reshaping(x_CG_Poisson,N,'CG solution');
%electric field

[Ex,Ey]=Efield(N,x_ref_resh);
xaxis=linspace(0,N,N); yaxis=xaxis;
figure
    quiver(xaxis,yaxis,Ex,Ey)
    xlabel('$x$ [a.u]', 'interpreter', 'latex')
    ylabel('$y$ [a.u]', 'interpreter', 'latex')
    title('Electric field', 'interpreter', 'latex')
    set (gca, 'fontsize', 22)
%Jacobi preconditioning
M=diag(diag(A));
[x_SD_Jacobi,countSD_Jacobi,errorSD_Jacobi]=solve_SD(M\A,M\rho);
[x_CG_Jacobi,countCG_Jacobi,errorCG_Jacobi]=solve_CG(M\A,M\rho);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Non linear conjugate gradient method
%cases N=2,3,4,5,6
clear,clc;
N=5; epsilon=1; sigma=1;
%x_init=randi([-10,10],[3*N,1])/10;
x_init=[0.5;0.3;0;-0.5;0.3;0;-0.5;-0.3;0;0.5;-0.3;0;0;0;0.6];
x_min=non_linCG(x_init);
Emin=energy(x_min);
%%
d=visualize_molecule(reshape(x_min,[3,N])',1.15,1,Emin);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




