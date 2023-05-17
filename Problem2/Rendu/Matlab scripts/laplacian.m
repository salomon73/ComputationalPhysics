N = 10; % Change this value to set the grid size

% Set up the grid
x = linspace(-1, 1, N);
[X,Y] = meshgrid(x);

% Set the grid spacing
dx = x(2) - x(1);

% Construct the Laplacian matrix
I = speye(N);
Dxx = spdiags([ones(N,1), -2*ones(N,1), ones(N,1)], [-1,0,1], N, N) / dx^2;
Dyy = spdiags([ones(N,1), -2*ones(N,1), ones(N,1)], [-1,0,1], N, N) / dx^2;
A = kron(Dxx,I) + kron(I,Dyy);

% Plot the sparsity pattern of the matrix
spy(A);

