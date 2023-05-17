N = 100; % Change this value to set the grid size
r = 0.5; % Change this value to set the circle radius

% Set up the grid
x = linspace(-1, 1, N);
[X,Y] = meshgrid(x);

% Find the points inside the circle
inside = (X.^2 + Y.^2) <= r^2;

% Plot the circle and the points inside it
figure;
hold on;
axis equal;
%viscircles([0,0], r, 'LineStyle', '-');
plot(X(inside), Y(inside), 'o');
xlim([-1 1])
ylim([-1 1])