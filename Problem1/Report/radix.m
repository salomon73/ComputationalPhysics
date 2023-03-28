function result=myfft(f)% f is an array input
    result = zeros(size(f));
    N = length(f);
    
    r=log2(N); % checking whether N is a power of two or not. Mandatory for
    % this algo.
    if abs(r-uint64(r))>1e-5 %arbitrary threshold 
        error('N has to be a power of 2.');
    end

    if N==1
        result = f;%this is computed directly from the formula
    elseif N==2
            result(1) = f(1) + f(2);%this is computed directly from the 
            % formula
            result(2) = f(1) - f(2);
    else 
        %recursive part, divide to conquer
        e = myfft(f(1:2:end)); %even part
        o = myfft(f(2:2:end)); %odd part
        for k=0:N/2-1
            W = exp(-2*1i*pi*k/N);
            result(k+1) = e(k+1) + W*o(k+1);
            result(k+1+N/2) = e(k+1)-W*o(k+1); %as usual, careful with the 
            % indices in Matlab, not the same as python
        end
    end
end