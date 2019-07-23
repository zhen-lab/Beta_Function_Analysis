function output = max_intensity_z_side_camera_single(vol)
% Create a maximum intensity projection in the z-direction (xy-plane)

% Take the maximum intensity projection.

im = squeeze(max(vol,[],3));


% Transpose and flip up/down
output = flipud(im);