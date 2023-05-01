function discrete_fourier = mydft(F)
    N    = length(F);   % spatial domain 
    temp = zeros(size(F));
    for m = 0:(N-1)     % Fourier space variable
        for n = 0:(N-1) % real space variable associated to domain of F
            temp(m+1)= temp(m+1) + F(n+1)*exp(-2*pi*1j*m*n/N);
        end
    end
    discrete_fourier = temp;
end