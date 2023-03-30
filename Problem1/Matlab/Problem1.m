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
transform = deltaT*(mydft((F)));
transform_mat = deltaT*(fft(F));
err = abs((transform-transform_mat)/transform_mat);

%% real and imaginary values of the Fourier transform of FID %%
figure
  subplot(2,1,1)
    plot(deltan, real(transform_mat), 'linewidth', 1);
    hold on
    h1 = plot(deltan, real(transform), 'linewidth', 1);
    xlim([min(deltan), max(deltan)])
    ylim([-30 300])
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
    ylim([-1e4, 8e4])
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
power = abs(transform).^2;
maxima = find(islocalmax(power, 'MinProminence', 1e4))
ind  = maxima(1)-100:maxima(1)+100;
ind2 = maxima(2)-100:maxima(2)+100;
grid = deltan(ind);

int_peak1 = trapz(grid,power(ind));
int_peak2 = trapz(deltan(ind2), power(ind2));

nprot1 = sqrt(int_peak1)
nprot2 = sqrt(int_peak2)
ratio = nprot1/nprot2


%% Test Mydft/Myfft

nrun = 10;
time = zeros(2,nrun);
for ii = 1:nrun
    tic 
        resdft = mydft(F);
    time(1,ii) = toc;
    tic 
        resfft = myfft(F);
    time(2,ii) = toc;
end
%%
figure
hold on 
plot(linspace(1,nrun, nrun),time(1,:), 'ro')
hold on 
plot(linspace(1,nrun, nrun),time(2,:), 'ko')
hold on
%plot(linspace(1,nrun,nrun), (log(N)/N)*ones(1,nrun), 'k-', 'linewidth', 1)
grid on
set(gca, 'fontsize', 22)
xlabel('$n_{run}$', 'interpreter', 'latex', 'fontsize',  28)
ylabel('$t$ [s]', 'interpreter', 'latex', 'fontsize', 28)
%% Test Mydft Myfft 2
nrun = 12;
time = zeros(2,nrun-4);
for ii = 5:nrun
    tic 
        resdft = mydft(F(1:2^ii));
    time(1,ii) = toc;
    tic 
        resfft = myfft(F(1:2^ii));
    time(2,ii) = toc;
end

%%
figure
hold on 
loglog(linspace(1,nrun, nrun),time(1,:), 'ro')
hold on 
loglog(linspace(1,nrun, nrun),time(2,:), 'ko')
hold on
%plot(linspace(1,nrun,nrun), (log(N)/N)*ones(1,nrun), 'k-', 'linewidth', 1)
grid on
set(gca, 'fontsize', 22)
xlabel('$n_{run}$', 'interpreter', 'latex', 'fontsize',  28)
ylabel('$t$ [s]', 'interpreter', 'latex', 'fontsize', 28)


%% Fourier Ptychography %%
addpath(genpath('Images'))
intens = imread('image_intensity.png');
phase  = imread('image_phase.png');
grey_int   = im2double(intens);
grey_phase = im2double(phase);

figure
    subplot(1,2,1)
      imshow(grey_int, [])
    subplot(1,2,2)
      imshow(grey_phase, [])

%% ptycho 2 %%

for ii=0%-2:1:2
    for jj=0%-2:1:2
        fname=['ptychography_' num2str(ii) '_' num2str(jj) '.png'];
        ptycho=fftshift(fft2(im2double(imread(fname))));
        ptycho_log=abs(ptycho); 
            imshow(log(ptycho_log),[]);
    end
end

%% ptycho 3 %%
%parameters
N      = 512; % Size of the image cutoff.png
guess  = ones(N);
fft2_guess   = fftshift(fft2(guess)); % fft on 2D array
cutoff = zeros(N);

%filter
for i=1:N
  for j=1:N
      if(i-(N-1)/2)^2+(j-(N-1)/2)^2 < 100^2
          cutoff(i,j)=1;
      else
          cutoff(i,j)=0;
      end
  end
end
filtred_guess=fft2_guess.*cutoff; % apply filter

%figure handling
xticklabel = [{}];
xticks = zeros(1,9);
for ii = 1:9
    xticklabel(ii) = num2cell(2^ii);
    xticks(ii)     = ii; 
end
yticklabel = xticklabel;

figure
  p1 = subplot(1,3,1);
    imagesc(log(abs(filtred_guess)))
    axis square
    colormap(gray);
    set(gca, 'fontsize', 22)
    xlabel('p', 'interpreter', 'latex', 'fontsize', 28)
    ylabel('p', 'interpreter', 'latex', 'fontsize', 28)

  p2 = subplot(1,3,2);
    imagesc(cutoff)
    colormap(gray);
    axis square
    set(gca, 'fontsize', 22)
    xlabel('p', 'interpreter', 'latex', 'fontsize', 28)
    ylabel('p', 'interpreter', 'latex', 'fontsize', 28)

  p3 = subplot(1,3,3);
    imagesc(log(abs(filtred_guess)))
    colormap(gray);
    axis square
    set(gca, 'fontsize', 22)
    xlabel('p', 'interpreter', 'latex', 'fontsize', 28)
    ylabel('p', 'interpreter', 'latex', 'fontsize', 28)

%% ptycho 4 %%
cutoff1 = im2double(imread('cutoff.png')); % filter read from cutoff.png
cutoff  = zeros(512,512);   % 'hand' built filter
for ii=1:N
  for jj=1:N
      if(ii-(N-1)/2)^2+(jj-(N-1)/2)^2 < 100^2
          cutoff(ii,jj)=1;
      else
          cutoff(ii,jj)=0;
      end
  end
end

guess   = fftshift(fft2(ones(512,512)));
image_filtered = guess.*cutoff;
figure
    imshow(log(image_filtered), []) % control everything is okay
inv = ifft2(image_filtered);
phase_inv = angle(inv);
phase_inv = sqrt(im2double(imread('ptychography_0_0.png'))).*exp(1j*phase_inv);
new_fft   = fft2(phase_inv);
new_fft   = new_fft.*cutoff;
figure
    imshow(log(new_fft), [])
figure
    imshow((ifft2(new_fft)), [])



%% Ptycho algorithm %%

guess = ones(512,512);
cutoff  = zeros(512,512);   % 'hand' built filter
for ii=1:N
  for jj=1:N
      if(ii-(N-1)/2)^2+(jj-(N-1)/2)^2 < 100^2
          cutoff(ii,jj)=1;
      else
          cutoff(ii,jj)=0;
      end
  end
end

%%
%Intens  = sqrt(im2double(imread('ptychography_0_0.png')));
Intens = sqrt(abs())

I_F     = fftshift(fft2(guess));
I_F_cut = I_F.*cutoff;
I_cut   = ifft2(I_F_cut);
phase   = angle(I_cut);
guess     = Intens.*exp(1j*phase);

new_F   = fftshift(fft2(new));
new_F   = new_F.*cutoff;

figure
    imshow(log(guess), [])

figure
   imshow(log(new_F), [])


%% ptycho 5 %% 
for ii =1:20
    for jj =1:21
        new_cutoff = circshift(cutoff,[50*ii, 50*(jj-1)]);
        im_filt = guess.* new_cutoff;
        inv = ifft2(im_filt);
        phase = angle(inv);
        phase = sqrt(im2double(imread('ptychography_0_0.png'))).*exp(1j*phase);
        new_fft   = fft2(phase);
        guess   = new_fft.*new_cutoff;
        %imshow(new_fft, [])
        %imshow(log(ifft2(new_fft)), [])
    end
end 
new_fft= guess;
figure
imshow(log((new_fft)), [])


%%
guess   = fftshift(fft2(ones(512,512)));
image_filtered = guess.*cutoff;


for ii = -2:1:2
    for jj = -2:1:2
        inv = ifft2(image_filtered);
        phase_inv = angle(inv);
        phase_inv = sqrt(im2double(imread(['ptychography_',num2str(ii),'_', num2str(jj),'.png']))).*exp(1j*phase_inv);
        new_fft   = fft2(phase_inv);
        new_fft   = log(new_fft.*cutoff);
   
        imshow(new_fft, [])
    end 
end







%% ptychography 4

inverseF_guess=ifft2(filtred_guess);
phase=angle(inverseF_guess);
magnitude=sqrt(im2double(imread('ptychography_0_0.png')));
ptycho=magnitude.*exp(1j*phase);
F_phase=fftshift(fft2(ptycho)); %FT of the complex created above
filtred_guess=F_phase;

%Applying the filter
F_phase_cutoff=F_phase.*cutoff;
figure
    imshow(log(F_phase_cutoff),[])

