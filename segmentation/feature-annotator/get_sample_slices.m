function y = get_sample_slices(x)
% y = GET_SAMPLE_SLICES(x)
%
%   This returns a list of values obtained by uniformly extracting slices
%   of the array x along the last dimension

s = size(x);

if s(end) <= 10
    
    y = x(:);
    
else
    
    step = floor(s(end)/10);
    idx = 1:step:s(end);
    
    y = [];
    
    for i = idx
        
        new_vals = get_slice(x, i);
        y = [y; new_vals(:)];
        
    end
    
end