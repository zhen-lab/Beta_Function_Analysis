function c_mask = circle_mask(ix, iy, cx, cy, r)

    % Width, height, center x, center y, radius
    [x, y] = meshgrid(-(cx - 1):(ix - cx), -(cy - 1):(iy - cy)); 

    % Points within the circle
    c_mask = ((x .^ 2 + y .^ 2) <= r ^ 2);

end