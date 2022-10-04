function [xyz_xy, pathname] = read_motorstage_position_uiget()

[filename, pathname]  = uigetfile({'*.tif'}, 'Select TIFF file');
xyz_col = read_motorstage_position(filename);
% xyz_xy = xyz_col(:, 1:2);

% xyz_col is a cell array
% Not sure if there is a better way to convert????
xyz_xy = zeros(size(xyz_col, 1), 2);
for j = 1:size(xyz_col,1)
    for i = 1:2
        xyz_xy(j,i) = str2double(xyz_col{j,i});
    end
end

end