function [ix, iy, valid] = pixelIndex(x, y, n_d, l_d)

    ix = floor(x/l_d + n_d/2) + 1;
    iy = floor(y/l_d + n_d/2) + 1;

    valid = ix >= 1 & ix <= n_d & iy >= 1 & iy <= n_d;

    ix(~valid) = NaN;
    iy(~valid) = NaN;
end