%% Gradient with finite differences %%
function g = grad(coord)
    N=length(coord);
    g=zeros(N,1); h=1e-3; id=eye(N);
    for i=1:N
        e=id(:,i);
        g(i)=(energy(coord+h*e)-energy(coord-h*e))/(2*h);
    end
end