%% myfft( F )
% =============================================
% 
% Computes the discrete Fourier transform of the
% function specified by the array F, following 
% the Radix2 algorithm.
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
% Written by S.Guinchard (03/23/23)   %
% ------------------------------------%
function discrete_fourier = myfft(F)

    discrete_fourier = zeros(size(F));
    N = length(f);
    r = log2(N); % checking whether N is a power of two or not. Mandatory for
                 % this algo.
    if abs(r-uint64(r))>1e-5 %arbitrary threshold 
        error('N has to be a power of 2.');
    end

    if N==1
            discrete_fourier = F; %this is computed directly from the formula
    elseif N==2
            discrete_fourier(1) = F(1) + F(2); %this is computed directly from the 
            discrete_fourier(2) = F(1) - F(2); % formula
    else 
        %recursive part
        e = myfft(F(1:2:end)); %even part
        o = myfft(F(2:2:end)); %odd part

        for k=0:N/2-1
            W = exp(-2*1i*pi*k/N);
            discrete_fourier(k+1) = e(k+1) + W*o(k+1);
            discrete_fourier(k+1+N/2) = e(k+1)-W*o(k+1);  
        end
    end
end