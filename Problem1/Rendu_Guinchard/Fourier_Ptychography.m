%% Fourier ptychography %%
% Salomon Guinchard      % 
% 03/30/2023             %
%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1 - Display initial intensity and phase %%
Intensity= im2double(imread("image_intensity.png"));
Phase    = im2double(imread("image_phase.png"));
figure
  subplot(1,2,1)
    imshow(Intensity,[]); 
  subplot(1,2,2)
    imshow(Phase,[]);

%% Microscope images %%

for ii=-2:1:2
    for jj=-2:1:2
        filename=['ptychography_' num2str(ii) '_' num2str(jj) '.png'];
        ptycho_F=fftshift(fft2(im2double(imread(filename))));
        ptycho_log=log(abs(ptycho_F));
        imshow(ptycho_log,[]);
    end
end

%% 2 - Ptychography_0_0.png %%

ptycho0 = fftshift(fft2(im2double(imread('ptychography_0_0.png'))));
figure
  imshow(log(ptycho0),[]);

%% 3 - Phase retrieval algorithm %% 

guess   = ones(size(Intensity));
guess_F = fftshift(fft2(guess)).*cutoff(0,0,1,1);
circle  = zeros(size(Intensity));
N       = length(Intensity(1,:));

% define circulartfilter
for i=1:N
  for j=1:N
      if(i-(N-1)/2)^2+(j-(N-1)/2)^2 < 100^2
          circle(i,j)=1;
      else
          cirlce(i,j)=0;
      end
  end
end

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
    imagesc(log(abs(guess_F)))
    axis square
    colormap(gray);
    set(gca, 'fontsize', 22)
    xlabel('p', 'interpreter', 'latex', 'fontsize', 28)
    ylabel('p', 'interpreter', 'latex', 'fontsize', 28)

  p2 = subplot(1,3,2); % Show the filter
    imagesc(circle)
    colormap(gray);
    axis square
    set(gca, 'fontsize', 22)
    xlabel('p', 'interpreter', 'latex', 'fontsize', 28)
    ylabel('p', 'interpreter', 'latex', 'fontsize', 28)

  p3 = subplot(1,3,3);
    imagesc(log(abs(guess_F)))
    colormap(gray);
    axis square
    set(gca, 'fontsize', 22)
    xlabel('p', 'interpreter', 'latex', 'fontsize', 28)
    ylabel('p', 'interpreter', 'latex', 'fontsize', 28)

%% 4 - Guess_function improvement %%

guess_F = fftshift(fft2(guess)); % take FFT of initial guess
iteration(0,0,512,guess_F);      % do one iteration

%% 5 - Repeat step 4 for all microscope images %%

j1 = -2:1:2; 
k1 = j1; 
new_guess=iteration(j1,k1,512,guess_F);
inverse(new_guess);

%% 6 - 20 Fourier ptychography iterations %%

for ii=1:20
    final_guess=iteration(j1,k1,512,guess_F);
    guess_F=final_guess;
end
inverse(final_guess);

%% Function definitions %%

function guess_updated = iteration(j1,k1,N,initial_guess)

    microscop    = cell(length(j1),length(k1));
    magnitude_   = cell(length(j1),length(k1));
    new_complex_ = cell(length(j1),length(k1));
    fft_improved = cell(length(j1),length(k1));
    negative_cutoff = cell(length(j1),length(k1));
    guess_updated   = zeros(N);

    for i=1:length(j1)
        for j=1:length(k1)
            fft_guess = initial_guess.*cutoff(j1,k1,i,j); 
            ifft_guess = ifft2(ifftshift(fft_guess));
            phase=angle(ifft_guess);
            filename = ['ptychography_' num2str(j1(i)) '_' num2str(k1(j)) '.png'];
            microscop{i,j}    = im2double(imread(filename));
            magnitude_{i,j}   = sqrt(microscop{i,j});
            new_complex_{i,j} = magnitude_{i,j}.*exp(1j.*phase);
            fft_improved{i,j} = fftshift(fft2(new_complex_{i,j})).*cutoff(j1,k1,i,j);
            negative_cutoff{i,j} = (1-cutoff(j1,k1,i,j)).*initial_guess;
            guess_updated = fft_improved{i,j}+negative_cutoff{i,j};
            initial_guess = guess_updated;
%             figure
%             imshow(log(abs(guess_updated)),[])
        end
    end
end

function [intensity,phase]=inverse(new_guess)

    ifft_guess = ifft2(ifftshift(new_guess));
    intensity  = abs(ifft_guess);
    phase      = angle(ifft_guess);
    figure
      subplot(1,2,1)
        imshow((intensity),[]);
      subplot(1,2,2)
        imshow(phase,[]);
end

function filter = cutoff(i1,j1,i,j)
    N = 512; 
    d = [50,50]; %displacement vector
    filter = zeros(N); 
    for m=1:N
        for n=1:N
            if (m-((N-1)/2+i1(i)*d(1)))^2+(n-((N-1)/2+j1(j)*d(2)))^2<100^2
                filter(m,n)=1;
            else
                filter(m,n)=0;
            end
        end
    end
end








