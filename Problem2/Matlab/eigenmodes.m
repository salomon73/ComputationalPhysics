%%       Eigenmodes of a vibrating string           %%
%   This script solves the eigenmodes of             %
%   a vibrating string, solving the wave equation    %
%   with Dirichlet boundary counditions              % 
% ____________________________________________________
%                                                    %
%        Written by S. Guinchard - 05/03/2023        %
%____________________________________________________%

% Input parameters
    l0 = 0;
    lN = 1;         % length of the rope 
    Npoints = 100;  % Number of meshpoints 
    dx = lN/(Npoints-1); 

% negative Laplace operator - matrix representation %%  
    A = zeros(Npoints, Npoints);
    for ii = 1:Npoints
        for jj = 1: Npoints
            if ii == jj 
                A(ii,jj) = 2/(dx*dx);
            elseif or((jj == ii+1) ,(jj == ii-1))
                A(ii,jj) = -1/(dx*dx);
            end
        end
    end
    if Npoints <=6 %shows result on console if not too big
        A
    end

% Amplitudes - lowest energy state
    w = rand(Npoints,1);
    w(1) = 0; w(end) = 0; % Impose Dirichlet boundary condition
    [eigenmode, eigenval] = eig_ipower(A,(pi/lN)^2);

% Amplitudes - 4 lowest energy states above ground
    for ii = 2:5
      [eigenv, eigenva] = eig_ipower(A,(ii*pi/lN)^2);
       eigenvals(ii-1)    = eigenva; 
       eigenvects(ii-1,:) = eigenv;
    end

%% Plot of results
    figure
        hold on 
        plot(linspace(0,1, Npoints),eigenmode, 'k-', 'linewidth', 1)
        set(gca, 'fontsize', 22)
        xlabel('x', 'interpreter', 'latex', 'fontsize', 28)
        ylabel('w(x)', 'interpreter', 'latex', 'fontsize', 28)
        legend(strcat('$\lambda=$',num2str(eigenval)), 'interpreter', 'latex', ...
            'location', 'northwest')

%% Multiplot for modes above ground state
    figure
        subplot(2,2,1)
            hold on 
            plot(linspace(0,1, Npoints),eigenvects(1,:), 'k-', 'linewidth', 1)
            set(gca, 'fontsize', 22)
            xlabel('x', 'interpreter', 'latex', 'fontsize', 28)
            ylabel('w(x)', 'interpreter', 'latex', 'fontsize', 28)
            legend(strcat('$\lambda=$',num2str(eigenvals(1))), 'interpreter', 'latex', ...
                'location', 'northwest')
        subplot(2,2,2)
            hold on 
            plot(linspace(0,1, Npoints),eigenvects(2,:), 'k-', 'linewidth', 1)
            set(gca, 'fontsize', 22)
            xlabel('x', 'interpreter', 'latex', 'fontsize', 28)
            ylabel('w(x)', 'interpreter', 'latex', 'fontsize', 28)
            legend(strcat('$\lambda=$',num2str(eigenvals(2))), 'interpreter', 'latex', ...
                'location', 'southwest')
        subplot(2,2,3)
            hold on 
            plot(linspace(0,1, Npoints),eigenvects(3,:), 'k-', 'linewidth', 1)
            set(gca, 'fontsize', 22)
            xlabel('x', 'interpreter', 'latex', 'fontsize', 28)
            ylabel('w(x)', 'interpreter', 'latex', 'fontsize', 28)
            legend(strcat('$\lambda=$',num2str(eigenvals(3))), 'interpreter', 'latex', ...
                'location', 'northeast')
        subplot(2,2,4)
            hold on 
            plot(linspace(0,1, Npoints),eigenvects(4,:), 'k-', 'linewidth', 1)
            set(gca, 'fontsize', 22)
            xlabel('x', 'interpreter', 'latex', 'fontsize', 28)
            ylabel('w(x)', 'interpreter', 'latex', 'fontsize', 28)
            legend(strcat('$\lambda=$',num2str(eigenvals(4))), 'interpreter', 'latex', ...
                'location', 'northeast')