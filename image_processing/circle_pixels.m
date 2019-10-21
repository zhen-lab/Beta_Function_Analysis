function [circle] = circle_pixels(width, height, b, a, r)

    [x, y] = meshgrid(1:height, 1:width);
    
    circle = ((x - a) .^ 2 +(y - b) .^ 2 <= r ^ 2);

end