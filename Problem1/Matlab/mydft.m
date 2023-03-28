%% mydft( F )
% =============================================
% 
% Computes the discrete Fourier transform of the
% function specified by the array F.
% 
% INPUT
% -----
%   -F:  size 1xlength(F) array containing the values
%        of the function F on its given domain 
%
% OUTPUT
% ------
%   -discrete_fourier: discrete Fourier transform
%                      of F
% ------------------------------------%
% Written by S.Guinchard (03/17/23)   %
% ------------------------------------%
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