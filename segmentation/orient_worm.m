function orient_worm(strain, animal)

% Set up directories
file_pattern = sprintf('%s_0%i (1).mat', strain, animal);
root = fullfile('A:', strain); %folder containing dataset
data_file = fullfile(root, file_pattern);

% Set a location for processed data
[~, f, ~] = fileparts(data_file);
    data_root = fullfile(root, f);
make_directory(data_root);
cd(root);

% Create a subfolder to hold useful summary data.
assets = fullfile(data_root, 'assets');
make_directory(assets);

% Create a subfolder to hold TIF files
red_stack_tiff_directory = fullfile(data_root, 'red');
green_stack_tiff_directory = fullfile(data_root, 'green');
make_directory(red_stack_tiff_directory);
make_directory(green_stack_tiff_directory);
D = red_stack_tiff_directory;

% acquire arraywz
mfile = matfile(data_file);
raw_buffer = mfile.images;
raw_buffer_trimmed = trim_movie(raw_buffer);

% This converts raw buffer images from the camera to correctly-oriented
% images.
get_zyla_vol = @(x) flip(x, 1);
get_top = @(x) x(1:256, :, :);
get_bottom = @(x) x(257:512, :, :);

% These two change depending on how the images are aligned on the camera,
% and need to be manually confirmed for every dataset (or at least every
% run).
get_green_slices = get_top;
get_red_slices = get_bottom;

raw_image = LazyArray(raw_buffer_trimmed, get_zyla_vol);
raw_green = LazyArray(raw_image, get_green_slices);
raw_red = LazyArray(raw_image, get_red_slices);

raw_green_z = LazyArray(raw_green, @max_intensity_z);
raw_red_z = LazyArray(raw_red, @max_intensity_z);

r0 = get_slice(raw_red_z, 1);
g0 = get_slice(raw_green_z, 1);

[optimizer, metric] = imregconfig('multimodal');
optimizer.InitialRadius = optimizer.InitialRadius/3.5; % Reduce the step size
tform = imregtform(r0, g0, 'translation', optimizer, metric);
rotate_red = @(x) imwarp(x, tform, 'OutputView', imref2d(size(r0)));

% r0_transformed = rotate_red(r0);

green = raw_green;
red = LazyArray(raw_red, rotate_red);

green_scaled = RangeScaledArray(green);
red_scaled = RangeScaledArray(red);

data = ColorArray({red_scaled, green_scaled}, 4, eye(2), false);
data_rgb = ColorArray({red_scaled, green_scaled}, 4, eye(3), true);
data_rgb_z = LazyArray(data_rgb, @max_intensity_z);


% Determine orientation
figure(3); clf;
imagesc(get_slice(data_rgb_z,1));
title('Select three points to determine orientation');
xlabel('(1) ventral nerve cord, (2) dorsal side, opposite first point, (3) anterior tip of head');

[x,y] = ginput(3);
tail_center = [mean(y(1:2)), mean(x(1:2))];
head_center = [y(3), x(3)];
center = mean([tail_center; head_center]);
angle = atan2(head_center(1)-tail_center(1), ...
    head_center(2) - tail_center(2))+pi;

%determine whether flipped
flipped = false;
cent_to_vent = [y(1), x(1)] - center;
cent_to_dors = [y(2), x(2)] - center;
ventral_angle = atan2(cent_to_vent(2), cent_to_vent(1));
dorsal_angle = atan2(cent_to_dors(2), cent_to_dors(1));

t = affine2d(eye(3));
t_rot = rotate_tform(t, angle, [128, 128]);

data_rotated = ImwarpedArray(data);
data_rotated_rgb = ImwarpedArray(data_rgb);

data_rotated_test = ImwarpedArray.update_tform(data_rotated_rgb, 1, t_rot);

if angdiff(ventral_angle, dorsal_angle) < 0
    flipped = true;
    data_rotated_test = LazyArray(data_rotated_test, @flip);
end
data_rotated_test_z = LazyArray(data_rotated_test, @max_intensity_z);

%% align to previous movie
target_stack = load_tiff_stack(D, 1);
target = max_intensity_z(target_stack);
source = get_slice(get_slice(data_rotated_test_z,1),1); 
[optimizer, metric] = imregconfig('multimodal');
% optimizer.InitialRadius = optimizer.InitialRadius/3.5; % Reduce the step size
tform = imregtform(source, target, 'rigid', optimizer, metric);
rotate_source = @(x) imwarp(x, tform, 'OutputView', imref2d(size(target)));
guess = rotate_source(source);
figure(2); 
% subplot(131); imshow(target); subplot(132); imshow(source); subplot(133); 
imshowpair(target, guess);
title('Alignment result');

%save transformation
save_location = fullfile(assets, 'transformation');
save(save_location, 't_rot', 'flipped');

