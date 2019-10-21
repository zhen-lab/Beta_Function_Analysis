function y = max_intensity_except(A, d)
% y = MAX_INTENSITY_EXCEPT(A, d)
%
%   Returns the maxmium intensity projection of A along the 'd' dimension.
%
% z = MAX_INTENSITY_EXCEPT(A, 3)
%
%   Maximum intensity projection in the x-, y-, and t-directions.
%
% Author: Vivek Venkatachalam (vivekv2@gmail.com)


N = ndims(A);
y = A;

for i = 1:N
    if i ~=d
        y = max(y, [], i);
    end
end

y = squeeze(y);
