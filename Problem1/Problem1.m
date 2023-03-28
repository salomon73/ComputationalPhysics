%% Load data for FID analysis %%
addpath(genpath('Matlab'))
F0 = load('Data/FID.dat');
F = complex(F0(:,1), F0(:,2));
N = length(F(:,1));
deltaT = 83.2e-6;
time   = deltaT*linspace(0,N-1, N); 
delta_nu = 1/(N*deltaT);
fc = 1/(2*deltaT);
ii = 0;
% while (-fc+ii*delta_nu < fc-delta_nu)
%     ii = ii+1;
%     fn(ii+1) = -fc+ii*delta_nu;
% end 
for ii = 1:N
    fn(ii) = (-(N/2)+ii-1)*(2/N*fc);
end 

%% FID decay plot %%
figure 
  subplot(2,1,1)
    h1 = plot(time, F0(:,1));
    xlim([0 0.339])
    xlabel('t [s]' , 'interpreter', 'latex')
    ylabel('$\mathcal{R}$(FID)' , 'interpreter', 'latex')
    set (gca ,'fontsize', 22)
    set(h1(1), 'Color', '#0072BD')
  subplot(2,1,2) 
    h2 = plot(time, F0(:,2));
    xlim([0 0.339])
    xlabel('t [s]' , 'interpreter', 'latex')
    ylabel('$\mathcal{I}$(FID)' , 'interpreter', 'latex')
    set (gca ,'fontsize', 22)
    set(h2(1), 'Color', '#D95319')

%% Definition of the discrete frequencies %%

f0 = (1067.93);
nu = 800.224e6;

deltan    = (fn-f0)/nu*1e6;
transform = (mydft((F)));
transform_mat = (fft(F));
err = abs((transform-transform_mat)/transform_mat);

%% real and imaginary values of the Fourier transform of FID %%
figure
  subplot(2,1,1)
    plot(deltan, real(transform_mat), 'linewidth', 1);
    hold on
    h1 = plot(deltan, real(transform), 'linewidth', 1);
    xlim([min(deltan), max(deltan)])
    grid on
    xlabel('$\delta_n$' , 'interpreter', 'latex')
    ylabel('$\mathcal{R}$(FID) [a.u.]' , 'interpreter', 'latex')
    set (gca ,'fontsize', 22)
    %set(h1(1), 'Color', '#0072BD')
    legend('$\mathcal{R}$(fft(FID))', '$\mathcal{R}$(mydft(FID))' ,'Location','northwest','Interpreter','latex');
  subplot(2,1,2) 
  grid on
    plot(deltan, imag(transform_mat), 'linewidth', 1);
    hold on 
    h2 = plot(deltan, imag(transform), 'linewidth',1);
    xlim([min(deltan), max(deltan)])
    grid on
    xlabel('$\delta_n$' , 'interpreter', 'latex')
    ylabel('$\mathcal{I}$(FID) [a.u.]' , 'interpreter', 'latex')
    set (gca ,'fontsize', 22)
    set(h2(1), 'Color', '#D95319')
    legend('$\mathcal{I}$(fft(FID))', '$\mathcal{I}$(mydft(FID))' ,'Location','northwest','Interpreter','latex');

%% Comparison mydft and matlab function %%
    figure
    subplot(2,1,1)
    h1 = plot(deltan, abs(transform).^2, 'linewidth', 1);
    grid on
    xlim([min(deltan), max(deltan)])
    ylim([-5e11, 12e12])
    set (gca ,'fontsize', 22)
    xlabel('$\delta_n$' , 'interpreter', 'latex', 'fontsize', 28)
    ylabel('$\mathcal{P}$ [a.u.]' , 'interpreter', 'latex', 'fontsize', 28)
    set(h1(1), 'Color', '#0072BD')
    legend('$\mathcal{P}$(mydft(FID))' ,'Location','northwest','Interpreter','latex');
  subplot(2,1,2) 
    h2 = plot(deltan, err);
    xlim([min(deltan), max(deltan)])
    ylim([-5e-15, 7e-14])
    grid on
    set (gca ,'fontsize', 22)
    xlabel('$\delta_n$' , 'interpreter', 'latex', 'fontsize', 28)
    ylabel('$\epsilon_r$' , 'interpreter', 'latex', 'fontsize', 28)
    legend('$|\frac{mydft-fft}{fft}|$' ,'Location','northwest','Interpreter','latex');
    set(h2(1), 'Color', '#D95319')
%%


