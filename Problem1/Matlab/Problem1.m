%% Load data for FID analysis %%
F0 = load('FID.dat');
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
    xlabel('t [s]' , 'interpreter', 'latex')
    ylabel('$\mathcal{R}$(FID)' , 'interpreter', 'latex')
    set (gca ,'fontsize', 22)
    set(h1(1), 'Color', '#0072BD')
  subplot(2,1,2) 
    h2 = plot(time, F0(:,2));
    xlabel('t [s]' , 'interpreter', 'latex')
    ylabel('$\mathcal{I}$(FID)' , 'interpreter', 'latex')
    set (gca ,'fontsize', 22)
    set(h2(1), 'Color', '#D95319')

%% Definition of the discrete frequencies %%

f0 = (1067.93);
nu = 800.224e6;

deltan    = (fn-f0)/nu*1e6;
transform = fftshift(mydft((F)));
transform_mat = fftshift(fft(F));

%% real and imaginary values of the Fourier transform of FID %%
figure
  subplot(2,1,1)
    h1 = plot(deltan, real(transform));
    xlabel('' , 'interpreter', 'latex')
    ylabel('$\mathcal{R}$(FID)' , 'interpreter', 'latex')
    set (gca ,'fontsize', 22)
    set(h1(1), 'Color', '#0072BD')
  subplot(2,1,2) 
    h2 = plot(deltan, imag(transform));
    xlabel('t [ms]' , 'interpreter', 'latex')
    ylabel('$\mathcal{I}$(FID)' , 'interpreter', 'latex')
    set (gca ,'fontsize', 22)
    set(h2(1), 'Color', '#D95319')

%% Comparison mydft and matlab function %%
    figure
    subplot(2,1,1)
    h1 = plot(deltan, abs(transform).^2);
    xlabel('$\delta$' , 'interpreter', 'latex')
    ylabel('P' , 'interpreter', 'latex')
    set (gca ,'fontsize', 22)
    set(h1(1), 'Color', '#0072BD')
  subplot(2,1,2) 
    h2 = plot(deltan, abs(transform_mat).^2);
    xlabel('$\delta$' , 'interpreter', 'latex')
    ylabel('P' , 'interpreter', 'latex')
    set (gca ,'fontsize', 22)
    set(h2(1), 'Color', '#D95319')
%%


