function raw_buffer_trimmed = trim_movie(raw_buffer)
% Finds first bright frame, and trims movie to duration of illumination

num_frames = 200; % number of frames for which the laser is on (e.g. 200)

% Finds peak in derivative of mean pixel value of each frame
mean_pix_val_array = mean(mean(mean(raw_buffer(:,:,:,:))));
mean_pix_val = squeeze(mean_pix_val_array(1,1,1,:));
[~, t] = max(diff(mean_pix_val));
first_bright_frame = t+1;

% Trims dataset to num_frames, beginning with first_bright_frame
bright_frames = first_bright_frame:2:first_bright_frame + num_frames - 1;
raw_buffer_trimmed = raw_buffer(:,:,:,bright_frames);