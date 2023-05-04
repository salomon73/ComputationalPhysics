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
    c = 1;    % ellipticity parameter





function H = hamiltonian(psi,pot)
% Apply the halmiltonian to the wavefunction psi
% with the potential V
me   = 9.10938291e-31;
hbar = 1.05457172647e-34;




end