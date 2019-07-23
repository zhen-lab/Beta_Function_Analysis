function [y, idx] = unpad(x, d)
% y = UNPAD(x)
%
%   Removes zeros along all dimensions.
%
% [y, idx] = UNPAD(x)
%
%   Returns an index cell array into the unpadded center of array x.
%
% y = UNPAD(x, d)
%
%   Unpads along the specified dimensions.
%
% Author: Vivek Venkatachalam (vivekv2@gmail.com)


N = ndims(x);
S = size(x);

if nargin < 2
    d = 1:N;
end

idx = {};
for i = d

    a = max_intensity_except(x, i);

    first = find(a, 1, 'first');
    last = find(a, 1, 'last');

    idx{i} = [first:last];

end

y = x(idx{:});