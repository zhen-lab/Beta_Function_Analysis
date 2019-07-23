function output = max_intensity_x_side_camera_single(vol)
% Create a maximum intensity projection in the z-direction (xy-plane)


output = flipud(squeeze(max(vol,[],2)));

