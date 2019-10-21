%% Set up directories
tic;

% USER: Set these to directory of dataset
root = pwd; %folder containing dataset 
file_pattern = ...
    'AVA_ATR_01 (1).mat'; 
%filename of dataset (e.g. AVA_01 (1).mat)

data_file = fullfile(root, file_pattern);
mfile = matfile(data_file);

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
data_filename = fullfile(data_root, 'mamut_data');
save_directory = fullfile(D, 'individual_neurons');
make_directory(save_directory);

% Initialize feature annotator
global feature_annotator_data;
feature_annotator_data.segmentation_options = struct();

% Acquire raw buffer
raw_buffer = MatfileArray(data_file, 'images');
size_raw_buffer = size(raw_buffer);

disp('Directory set up.');


% Trim movie and select only odd numbers (1, 3, 5, ..., 197, 199)
% Finds first bright frame, and trims movie to duration of illumination
% Only odd numbers are used because volumes expand/shrink every other
% volume. Odd numbers or even numbers are aligned.

num_frames = 100; % number of frames for which the laser is on (e.g. 200), and only take odd numbers

% Finds peak in derivative of mean pixel value of each frame
mean_pix_val_array = mean(mean(mean(raw_buffer(:,:,:,:))));
mean_pix_val = squeeze(mean_pix_val_array(1,1,1,:));
[~, t] = max(diff(mean_pix_val));
first_bright_frame = t+1;

% Trims dataset to num_frames, beginning with first_bright_frame
bright_frames = first_bright_frame:2:first_bright_frame + num_frames*2 - 1;
% valid_frames = 1:2:200;
raw_buffer_cropped = TimeMappedArray(raw_buffer, bright_frames);

disp('Time window trimmed.');


%% Create an aligned array

% Reference time for global rotation alignment.
t0 = 1; % Frame in the middle appears to most stable for pattern.

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

raw_image = LazyArray(raw_buffer_cropped, get_zyla_vol);
raw_green = LazyArray(raw_image, get_green_slices);
raw_red = LazyArray(raw_image, get_red_slices);

raw_green_z = LazyArray(raw_green, @max_intensity_z);
raw_red_z = LazyArray(raw_red, @max_intensity_z);

r0 = get_slice(raw_red_z, t0);
g0 = get_slice(raw_green_z, t0);
figure(2); clf; subplot(121); imshowpair(g0, r0);

[optimizer, metric] = imregconfig('multimodal'); % 'monomodal' does not work well
optimizer.InitialRadius = optimizer.InitialRadius/3.5; % Reduce the InitialRadius to avoid non-converging issue
tform = imregtform(g0, r0, 'translation', optimizer, metric); % Move GFP
rotate_red = @(x) imwarp(x, tform, 'OutputView', imref2d(size(g0)));

g0_transformed = rotate_red(g0);
figure(2); subplot(122); imshowpair(g0_transformed, r0);

registration_filename = fullfile(assets, 'red_green_alignment.png');
saveas(gcf, registration_filename);

% green = raw_green;
% red = LazyArray(raw_red, rotate_red);
green = LazyArray(raw_green, rotate_red);
red = raw_red;

% green_scaled = RangeScaledArray(green, 3);
% red_scaled = RangeScaledArray(red, 2);
green_scaled = RangeScaledArray(green, 1);
red_scaled = RangeScaledArray(red, 1);

disp('Registration completed.');

%% Generate max intensity projections (MIPs)

pixel_ratio = [1, 1, 2];
red_mips = ThreeViewArray(red_scaled, pixel_ratio);
green_mips = ThreeViewArray(green_scaled, pixel_ratio);

red_mips_colored = ColorArray(red_mips, 3, [1;0;0], true);
green_mips_colored = ColorArray(green_mips, 3, [0;1;0], true);

mips_side_by_side = SideBySideArray(...
    {red_mips_colored, green_mips_colored}, 2, 10);
mips_side_by_side = AnnotatedArray(mips_side_by_side);
mips_side_by_side.Functions{1} = @AnnotatedArray.timestamp;
mips_side_by_side.FunctionArgs{1} = {[10, 10]};

%check projections
figure(3); clf; imshow(get_slice(mips_side_by_side, 10));

mips = ColorArray({red_mips, green_mips}, 3, eye(3), true);
figure(4); clf; imshow(autoscale(get_rgb_slice(mips, 1)));
figure(5); clf; imshow(autoscale(get_rgb_slice(mips, size(mips, 4))));

mips = ColorArray({red_mips, green_mips}, 3, [1,0; 0,1; 0,0], true);
mips = AnnotatedArray(mips);
mips.Functions{1} = @AnnotatedArray.timestamp;
mips.FunctionArgs{1} = {[10, 10]};

% Write MIPs to video
mipsfile_1 = fullfile(assets, 'MIP_side_by_side.mp4');
video_from_array(mips_side_by_side, mipsfile_1, ...
    'framerate', 10);

mipsfile_2 = fullfile(assets, 'MIP.mp4');
video_from_array(mips, mipsfile_2, ...
    'framerate', 10);

disp('MIP videos created.');


%% Determine orientation
% Select three points to configure orientation
% (1) ventral nerve cord
% (2) dorsal side directly opposite first point
% (3) anterior tip of head

data = ColorArray({red_scaled, green_scaled}, 4, eye(2), false);
data_rgb = ColorArray({red_scaled, green_scaled}, 4, eye(3), true);
data_rgb_z = LazyArray(data_rgb, @max_intensity_z);

figure(3); clf;
% imagesc(get_slice(data_rgb_z, 1)); 
imagesc(51+0.8*get_slice(data_rgb_z, 1)); % To see the cursor position
title('Select three points to determine orientation');
xlabel('(1) ventral nerve cord, (2) dorsal side, opposite first point, (3) anterior tip of head');

hold on;
num_pts = 3; x = zeros(num_pts, 1); y = zeros(num_pts, 1);
for np = 1:num_pts
    [x(np),y(np)] = ginput(1);
    plot(x(np), y(np), 'om'); text(x(np)+2, y(np)+2, num2str(np), 'color', 'm');
end
drawnow;
hold off;
% [x,y] = ginput(3);
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

if ventral_angle > dorsal_angle
    data_rotated_test = LazyArray(data_rotated_test, @flip);
end

data_rotated_test_z = LazyArray(data_rotated_test, @max_intensity_z);
figure(5); clf; 
imagesc(get_slice(data_rotated_test_z, 1)); 
title('Sample frame after rotation');
drawnow;
disp('If it looks good, apply transformation to all frames.');


%% If it looks good, apply transformation to all frames
% otherwise, repeat previous cell

for i = 1:num_frames
    data_rotated = ImwarpedArray.update_tform(data_rotated, i, t_rot);
    data_rotated_rgb = ImwarpedArray.update_tform(data_rotated_rgb, ...
        i, t_rot);
end

if ventral_angle > dorsal_angle
    data_rotated = LazyArray(data_rotated, @flip);
    data_rotated_rgb = LazyArray(data_rotated_rgb, @flip);
end

disp('Rotation applied to all frames.');


%% Crop movie to region of interest
% place one point above the animal, and one point below

figure(4); clf; 
% imagesc(get_slice(data_rotated_test_z, 1));
imagesc(51+0.8*get_slice(data_rotated_test_z, 1)); % To see the cursor position
title('Select two points to determine crop');
xlabel('(1) above highest neuron, (2) below lowerst neuron');

% [~,y] = ginput(2);
hold on;
num_pts = 2; y = zeros(num_pts,1);
for np = 1:num_pts
    [~,y(np)] = ginput(1);
    plot([1 size(g0,1)], [y(np) y(np)], 'm');
end
drawnow;
hold off;
bound_1 = round(y(1));
bound_2 = round(y(2));
if bound_1 < bound_2
    crop_frames = bound_1:bound_2;
else
    crop_frames = bound_2:bound_1;
end

data_trimmed = LazyArray(data_rotated, @(x) x(crop_frames, :, :, :, :));
data_trimmed_rgb = LazyArray(data_rotated_rgb, ...
    @(x) x(crop_frames, :, :, :, :));

data_trimmed_rgb_z = LazyArray(data_trimmed_rgb, @max_intensity_z);

figure(1); image(max_intensity_z(get_slice(data_trimmed_rgb, num_frames)));

for i = 1:10:size(data_trimmed_rgb_z, 4)
    imagesc(get_slice(data_trimmed_rgb_z, i));
    pause(0.01);
    drawnow;
end

disp('Movie trimmed.');

%% Filter and normalize

bump_finder = bump_detector([1,1], [3,3], [51 51]);
normalization = 0.7 * sum_all(bump_finder(bump_finder>0)*255);

dc = @(x) deconvlucy(x, fspecial('gaussian', 5, 2), 1);
flt = @(x, h) imfilter(double(dc(x)), h, 'circular');
flt_and_normalize = @(x) uint8(flt(x, bump_finder)/normalization*256-5);

flt_normalize_slice = @(x) flt_and_normalize(x(crop_frames,...
    :,:,:,:));

data_trimmed_filtered_normalized = LazyArray(data_trimmed, ...
    flt_and_normalize);
data_trimmed_filtered_normalized_rgb = LazyArray(data_trimmed_rgb, ...
    flt_and_normalize);

figure(4); imshow(max_intensity_z(get_slice(data_trimmed_rgb, num_frames)));
title('pre filt/deconv')
figure(5); imshow(max_intensity_z(get_slice(...
    data_trimmed_filtered_normalized_rgb, num_frames)));
title('post filt/deconv')

disp('Filtering completed.');
toc;

% Write TIFs
red_stack = LazyArray(data_trimmed_filtered_normalized, @(x) squeeze(x(:,:,:,1)));
green_stack = LazyArray(data_trimmed_filtered_normalized, @(x) squeeze(x(:,:,:,2)));

array_to_TIFF_stacks(red_stack, red_stack_tiff_directory);
array_to_TIFF_stacks(green_stack, green_stack_tiff_directory);


%%%%%%%%% Save to big data viewer array

data_for_mamut = LazyArray(data_trimmed_filtered_normalized,...
    @(x) permute(x(:,:,:,:,:), [2,1,3,4,5]));
write_big_data_viewer_array(data_filename, data_for_mamut);

disp('Tiff files saved, Mamut array created');
toc;

%% Automatic annotation

tic;
Dirc = dir([D, '\*.tif']);
num_frames = length(Dirc(not([Dirc.isdir])));
v = load_tiff_stack(D, 1);
fprintf('interpolating volume data \n');
% vq = v;
vq = interp3(double(v)); % interpolation of data to avoid consecutive maxima

fprintf('finding local maxima \n');
num_neurons_to_identify = 20;
maxima = find_local_maxima(v, [5, 5, 3], 'max', num_neurons_to_identify, ...
    'minimum_separation', [0.5, 0.5, 0.5]);
% maxima = find_local_maxima(vq, [5, 5, 3]);
maxima = ceil(maxima/2);
flength = size(maxima, 1);
fsize = [50, 50, 5];

features = cell(1, flength);
for i = 1:flength
    name = sprintf('neuron_%03d', i);
    f = new_feature(name, num_frames, fsize);
%     f = translate_feature(f, maxima(i,:));
    f = translate_feature(f, maxima(i,:) - fsize/2);
    f.modified_coordinates(1,:) = f.coordinates(1,:);
    features{i} = f;
end

fprintf('saving data \n');
feature_annotator_data.features = features;
save_features(features, D);
save_tracks = tracks_from_features(features);
mamut_xml_from_tracks(save_tracks, data_filename, [1 1 2]);
fprintf('automatic annotation completed. \n');
toc;

%% NOW ANNOTATE IN MaMuT

%% Import tracks from MaMuT

tracks = tracks_from_mamut(data_filename, [1 1 2]);
features = features_from_tracks(tracks, num_frames);

feature_annotator_data.features = features;
save_features(features, D);


%% Prepare for tracking

v1 = load_tiff_stack(D, 1);
N_features = length(features);
cutoff = 0.1;
max_tries = 30;
discount = 0.95;
t_0 = 1;
track_size = [25 25 5];

features = set_feature_size(features, track_size);

options = struct(...
    'drift', [2, 2, 1], ...
    'imreg_padding', [10, 10, 2]...
);

feature_annotator_data.features = features;
feature_annotator_data.segmentation_options = ...
    merge_struct(options, feature_annotator_data.segmentation_options);
feature_annotator_data.image_location = D;


%% Track neurons

tic;

parfor i = 1:length(features)
    t_ref_best = 1;
    f = features{i};
    disp(['Registering feature ', num2str(i)]);

    good_frames = get_modified_frames(f);

    qualities = zeros(1, num_frames);
    qualities(good_frames) = 1;
    
    depths = ones(1, num_frames)*1000;
    depths(good_frames) = 0;

    for t = 1:num_frames
        q_remaining = qualities .* discount.^(depths);
        t_remaining = 1:num_frames;
        tries = 0;
        modified_max_tries = min([t-1, max_tries]);
        q = 0;
        f0 = f; % feature(i)
        while q < cutoff && ~isempty(q_remaining) && ...
              tries < modified_max_tries

            idx_next = find(q_remaining == max(q_remaining), 1);
            
            t_ref = t_remaining(idx_next);
            
            t_remaining(idx_next) = [];
            q_remaining(idx_next) = [];
            
            f_new = register_feature_frame(f0, D, t, ... 
                't_ref', t_ref, options);
            q_new = f_new.registration_info{t}.image_registration_score;
            
            if q_new > q
                f = f_new;
                q = q_new;
                t_ref_best = t_ref;
            end

            tries = tries + 1;
            
        end

        fprintf('Time %d okay in %d tries. Score: %d. Parent: %d\n', ...
            t, tries, q, t_ref_best);
        
        depths(t) = depths(t_ref_best) + 1;        
        qualities(t) = q;

    end
    
    f.registration_summary.qualities = qualities;
    f.registration_summary.depths = depths;
    f.registration_summary.time_complete = now;

    features{i} = f;

    filename = fullfile(save_directory, sprintf('neuron_%04d.mat', i));
    save_features(f, filename);

end

feature_annotator_data.features = features;
save_features(features, D);


% Translate back to MaMuT

save_tracks = tracks_from_features(features);
mamut_xml_from_tracks(save_tracks, data_filename, [1 1 2]);

toc;


%% NOW PROOFREAD IN MaMuT


%% Re-import proofread tracks from MaMuT

tracks = tracks_from_mamut(data_filename, [1 1 2]);
features = features_from_tracks(tracks, num_frames);

feature_annotator_data.features = features;
save_features(features, D);


%% Extract calcium traces
% [Optional] reads intensities of each channel at each time point, and
% saves them as time series
tic;
feature_size = [3, 3, 4];
options = struct(...
    'npoints', 10 ...
);

% features = load_features(D);
features = set_feature_size(features, feature_size);

rfp = get_intensity(features, red_stack_tiff_directory, options);
gcamp = get_intensity(features, green_stack_tiff_directory, options);


filename = fullfile(data_root, 'traces.mat');
save(filename, 'gcamp', 'rfp');

figure(1); clf; plot(gcamp);
saveas(gcf, fullfile(data_root, 'gcamp.jpg'));
figure(2); clf; imagesc(transpose(gcamp));
saveas(gcf, fullfile(data_root, 'heatmap.jpg'));
toc;


%% Load features from files

features = load_features(D);

for i = 1:length(features)
    filename = fullfile(save_directory, sprintf('neuron_%04d.mat', i));
    if exist(filename, 'file')
        n = load_features(filename);
        features{i} = n{1};
    end
end