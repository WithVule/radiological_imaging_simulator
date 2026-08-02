function [x_Q, y_Q] = positionOnFocalPoint(xi, dim)

x_Q = (xi(:, 1) - 0.5).*(dim(1));
y_Q = (xi(:, 2) - 0.5).*(dim(2));
end