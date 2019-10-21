function [y, idx] = trim_edges(x, varargin)
% y = TRIM_EDGES(x)
%
%   Trims zeroes from the boundary of an ND array.
%
% [y, indx] = TRIM_EDGES(x)
%
%   Additionally returns a cell array of indices describing the trim.
%
% Author: Vivek Venkatachalam (vivekv2@gmail.com)

default_options = struct(...
    'filter', [] ...
);

input_options = varargin2struct(varargin{:}); 
options = mergestruct(default_options, input_options);

if ~isempty(options.filter)
    
    if isnumeric(options.filter)
        x_0 = imfilter(x, options.filter);
    else
        x_0 = options.filter(x);
    end
    
else
    
    x_0 = x;
    
end

[~, idx] = unpad(x_0);

y = x(idx{:});