%%           2D Elliptical quantum well             %%
%  This scripts  implements  the  solution of the    %
%  time  independent  Schrödinger  equation on a     %
%  two-dimensional cartesian grid, for a constant    %
%                  non zero potential                %
% ____________________________________________________
%                                                    %
%        Written by S. Guinchard - 05/03/2023        %
%____________________________________________________%


% Input parameters
    r = 1e-9; % r = 1nm
    a = 5*r;  % range over which the grid is defined
    N = 100;  % number of points in each dir 
    c  = 1;    % ellipticity parameter
    V0 = -1.5; % eV - potential well
    constants = [r, a, N, c];

% define the potential matrix
% It will be converted to sparse
% inside hamiltonian function
    pot = zeros(N,N);
    grid = linspace(-a/2, a/2,N);
    [X,Y] = meshgrid(grid);
    for ii = 1:N
        for jj = 1:N
            if (X(ii,jj)^2+Y(ii,jj)^2/c^2)<=r^2
                pot(ii,jj) = V0;
            end
        end 
    end 


% call Hamiltonian function to visualize 
% the halmiltonian
    Hm = hamiltonian_matrix(pot, constants);

% Need to compute the eigenvalues to find the bound states
    eigs(Hm)
    threshold = 0; % we want all negative eigenvalues
    
    lowest_energy = eigs(Hm,1);
    dE = abs(lowest_energy)/20;
    neg_eig = lowest_energy;

% Brute force, by hand
% They're all here! 
% Procedure: eigs(Hm,1,guessed_eigenval)
    negative_eigenvalues = [-49.6 -6.87 -4.05, -3.26 -2.47 -1.90 -1.54];


 %% Bound states for three lowest energies: 
[bound,val] = eigs(Hm,1,-6.87); % remplacer par les eigenvalues
                               % en revanche, pas sur de la conversion avec
                               % reshape
bound = reshape(bound, 100, 100);
figure
pcolor(ctranspose(bound)*bound) % plot densité de proba
shading interp

% Don't bother with that loop
% loop to find all the eigenvalues
% please help 
% 
%     for ii = 1:20
%         val = eigs(Hm,1,lowest_energy+ii*dE);
%         if abs(val-neg_eig(end))<1e-2
%         else 
%            cat(2,neg_eig,val);
%         end 
%     end


% function definitions below %
%____________________________%

% Hamiltonian operator matrix as a sparse
function Hm = hamiltonian_matrix(pot, constants)
    
    % input constants
    me   = 9.10938291e-31;
    hbar = 1.05457172647e-34;
    N = constants(3); 
    a = constants(2);


    % grid parameters
    x = linspace(-a/2, a/2, N);
    dx = x(2) - x(1);
    
    % Construct the Laplacian matrix - centered finite differences
    I = speye(N);
    Dxx = spdiags([ones(N,1), -2*ones(N,1), ones(N,1)], [-1,0,1], N, N) / dx^2;
    Dyy = spdiags([ones(N,1), -2*ones(N,1), ones(N,1)], [-1,0,1], N, N) / dx^2;
    A  = kron(Dxx,I) + kron(I,Dyy);
    Hm = -(hbar^2/(2*me))*A + kron(sparse(pot), I); 
    
    % Plot the sparsity pattern of the matrix
    % to visualize the Hamiltonian
    spy(Hm);


end

% This function needs a wave function defined on the grid
% to apply the hamiltonian
function H = hamiltonian_grid(psi,pot, constants)
    % Apply the halmiltonian to the wavefunction psi
    % with the potential V
    % The operator can be rewritten as a matrix, see function 
    % hamiltonian_matrix
    me   = 9.10938291e-31;
    hbar = 1.05457172647e-34;
    a = constants(2);
    N = constants(3);
    dx = a/(N-1);
    dy = dx; % square grid
    
    H = zeros(size(psi));
    % Apply hamiltonian operator inside the grid
    for ii = 2:N-1
        for jj = 2:N-1
            H(ii,jj) = -(hbar^2/2*me)*( -2*((dx^2 + dy^2)/(dx^2*dy^2))*psi(ii,jj) ...
                       + (1/dx)^2*(psi(ii-1,jj)+ psi(ii+1,jj)) ...
                       + (1/dy)^2*(psi(ii,jj-1)+ psi(ii,jj+1)) ...
                       + pot(ii,jj)*psi(ii,jj));
        end
    end
    % Apply hamiltonian operator on the boundary (forward) left and right 
    for jj = [1,N]
        for ii = 1:N
            H(ii,jj) = -(hbar^2/2*me)*(...
                        (psi(ii+2,jj) - 2*psi(ii+1, jj) + psi(ii,jj)) / (dx^2) ...
                      + (psi(ii,jj+2) - 2*psi(ii, jj+1) + psi(ii,jj)) / (dy^2) ...
                                     ) + pot(ii,jj)* psi(ii,jj);
        end
    end
    % Apply hamiltonian operator on the boundary (forward) up and down 
    for ii = [1,N]
        for jj = 2:N-1 % corners already done
            H(ii,jj) = -(hbar^2/2*me)*(...
                        (psi(ii+2,jj) - 2*psi(ii+1, jj) + psi(ii,jj)) / (dx^2) ...
                      + (psi(ii,jj+2) - 2*psi(ii, jj+1) + psi(ii,jj)) / (dy^2) ...
                                     ) + pot(ii,jj)* psi(ii,jj);
        end
    end
end

