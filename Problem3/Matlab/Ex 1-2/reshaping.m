%% Reshaping and ploting function %%
function X_reshaped = reshaping(X,N,name)
    X_reshaped=reshape(X,N,N);
    xaxis=linspace(0,N,N); yaxis=xaxis;
    figure
    p=pcolor(xaxis,yaxis,X_reshaped);
    shading interp
    colorbar
    p.FaceColor='interp';
    xlabel('$x$ [a.u]')
    ylabel('$y$ [a.u]')
    title(name)
    %contourf(X_reshaped); colorbar
end