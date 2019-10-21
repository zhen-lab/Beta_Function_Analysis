function output = max_intensity_z_top_camera(vol)
% Create a maximum intensity projection in the z-direction (xy-plane)

% Take the maximum intensity projection.

im = squeeze(max(vol(:,:,1,:),[],4));

output = im;