function output = max_intensity_x_top_camera(vol)
% Create a maximum intensity projection in the x-direction (yz-plane)

output = squeeze(max(vol(:,:,1,:),[],2));

