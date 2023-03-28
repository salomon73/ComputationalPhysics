fs = 100;               % sampling frequency
t = 0:(1/fs):(10-1/fs); % time vector
S = cos(2*pi*15*t);
n = length(S);
X = fft(S);
f = (0:n-1)*(fs/n);     %frequency range
power = abs(X).^2/n;    %power
plot(f,power)
Y = fftshift(X);
fshift = (-n/2:n/2-1)*(fs/n); % zero-centered frequency range
powershift = abs(Y).^2/n;     % zero-centered power
figure
plot(fshift,powershift)

%%
t=linspace(-4.0,4.0,N+1); 
f=gaussmf(t,[1 0]); 
g=fft(f);
plot(t,g);