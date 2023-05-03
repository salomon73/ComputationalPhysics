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
    lN = 1;      % length of the rope 
    Npoints = 100; % Number of meshpoints 

% negative Laplace operator - matrix representation %%  
    A = zeros(Npoints, Npoints);
    for ii = 1:Npoints
        for jj = 1: Npoints
            if ii == jj 
                A(ii,jj) = 1/(Npoints-1)^2*2;
            elseif or((jj == ii+1) ,(jj == ii-1))
                A(ii,jj) = -1/(Npoints-1)^2;
            end
        end
    end
    if Npoints <=6 %shows result on console if not too big
        A
    end

% Amplitudes - lowest energy state
    w = rand(Npoints,1);
    w(1) = 0; w(end) = 0; % Impose Dirichlet boundary condition
    [eigenmode, eigenval] = eig_rq(A,);




% 
%     
% % Amplitudes - above lowest (4)
% vector = [(pi/9)^5,(pi/9)^4.5,(pi/9)^4,(pi/9)^3.5];
% for ii = 1:4
%     w = rand(Npoints,1);
%     w(1) = 0; w(end) = 0; % Impose Dirichlet boundary condition
%     [eigenmode_, eigenval_] = eig_rq(A,(pi/9)^3.5);
% end
% Plot of results
    figure
    hold on 
    plot(linspace(0,1, Npoints),eigenmode, 'k-', 'linewidth', 2)
    set(gca, 'fontsize', 22)
    xlabel('x', 'interpreter', 'latex', 'fontsize', 28)
    ylabel('w(x)', 'interpreter', 'latex', 'fontsize', 28)

