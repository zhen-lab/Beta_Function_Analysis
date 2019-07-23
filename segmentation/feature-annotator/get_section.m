function section = get_section(image, start_indices, feature_size)
% function section = GET_SECTION(image, start_indices, feature_size)
%
%   Extracts a section of an image with coordinates starting at
%   'start_indices' and a size of 'size'.  Boundaries are padded with the 
%   average value on the face to allow extraction of regions that jut out
%   of the original image.
%
%   If feature_size has fewer dimensions than the image, no sectioning is
%   done on tailing dimensions.
%
% Author: Vivek Venkatachalam (vivekv2@gmail.com)

size_image = size(image);
N = length(size_image); 

M = length(start_indices);
assert(M==length(feature_size), ...
    'Start index and feature size must be the same dimension.');        
start_indices(M+1:N) = 0.5;
feature_size(M+1:N) = size_image(M+1:N);

% Determine the average value of border pixels.
border_pixels = [];
for i=1:N
    for j=1:N
        if j==i
            idx{j} = [1 size_image(j)];
        else
            idx{j} = [1:size_image(j)];
        end
    end
    border_pixels = [border_pixels ...
                     reshape(image(idx{:}),1,numel(image(idx{:})))];
end
background_intensity = mean(border_pixels);
clear border_pixels;

% Create a larger version of vol to allow features near edges. pad with
% the average value on the border pixels.

bigimage_start_indices = feature_size;
bigimage_end_indices = bigimage_start_indices + size_image - 1;

for i=1:N
    idx{i} = bigimage_start_indices(i):bigimage_end_indices(i);
end

image_full_idx = cell(1,N);
[image_full_idx{:}] = deal(':');

image_large_size = size(image) + 2*(feature_size-1);

image_large = ones(image_large_size, element_class(image)) ...
    * background_intensity;
image_large(idx{:})=image(image_full_idx{:});
clear idx image; % don't need this anymore

new_start_indices = start_indices + bigimage_start_indices - 1;
end_indices = new_start_indices + feature_size - 1;

for i=1:N
    idx{i} = round(new_start_indices(i):end_indices(i));
end
section = image_large(idx{:});
