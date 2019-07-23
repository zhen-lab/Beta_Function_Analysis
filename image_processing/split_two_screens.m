function [ left_cell, right_cell ] = split_two_screens( imagelist )

framenum = size(imagelist, 1);
width = size(imagelist{1,1}, 2);

left_cell = cell(framenum, 1);
right_cell = cell(framenum, 1);

for j = 1:framenum
    
    left = imagelist{j,1}(:, 1:width/2); % The left screen, RFP channel
    left_cell{j} = left;
    
    right = imagelist{j,1}(:, (width/2+1):end); % The right screen, GFP channel
    right_cell{j} = right;

end

end