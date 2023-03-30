function result=FFT(input_a)

 N = length(input_a);

 if N == 1
 result = input_a;
 else
 input_even = input_a(1:2:N); %even elements
 input_odd = input_a(2:2:N); %odd elements

 reseven = myfft(input_even);
 resodd = myfft(input_odd);

 result = input_a; %create a vector result, same size as input

 for k=1:N/2
 W = exp(-2*pi*1i*(k-1)/N);
 result(k) = reseven(k) + W*resodd(k);
 result(k+N/2) = reseven(k) - W*resodd(k);
 end

 end

 end