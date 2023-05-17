%% Concergence study %% 
function [Nmax,pi_num]=convergence(min,max)
N=min:2:max; Req=zeros(1,length(N));
diff=zeros(1,length(N));
for i=1:length(N)
    [A,b]=generate_ab(N(i));
    Proj=projection(N(i));
    projA=Proj*A*Proj';
    y=projA\(Proj*b); x=Proj'*y;
    [V,~,~]=reshaping(N(i),x);
    Req(i)=V((N(i)-1)/2,(N(i)-1)/2)-V((N(i)-1)/2+1,(N(i)-1)/2+2);
end
figure
plot(N,Req,'.');
xlabel('$N$')
ylabel('$R_{eq}$')
pi_num=(Req(end)+0.5)^(-1)*4;
for i=1:length(N)-1
    diff(i)=abs(Req(i+1)-Req(i));
    if abs(diff(i))<=1e-2
       Nmax=N(i+1);
       return;
    end
end