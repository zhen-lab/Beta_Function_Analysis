function y = get_full_path(x)
% y = GET_FULL_PATH(x)
%
%   Returns the full path of a (possibly relative) path x.
%
% Author: Vivek Venkatachalam (vivekv2@gmail.com)

current_directory = pwd;

if strfind(x, current_directory) % x is an absolute path
    
    y = x;
    
else % x is a relative path
    
    y = fullfile(current_directory, x);
    
end