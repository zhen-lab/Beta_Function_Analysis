function output = max_intensity_y_top_camera(vol)
% Create a maximum intensity projection in the y-direction (xz-plane)


output = squeeze(max(vol(:,:,1,:),[],1))';