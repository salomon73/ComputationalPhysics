%% Hessian matrix projected along a preferential direction d %%
function H = myhess(x,d)
    h=1e-3; d=d/norm(d); %normalize d
    H=(energy(x+h*d)-2*energy(x)+energy(x-h*d))/h^2;
end