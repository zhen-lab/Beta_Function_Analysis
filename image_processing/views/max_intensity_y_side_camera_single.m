function output = max_intensity_y_side_camera_single(vol)
% Create a maximum intensity projection in the z-direction (xy-plane)


output = squeeze(max(vol,[],1))';
