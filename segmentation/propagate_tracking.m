function propagate_tracking(strain, animal)

%% Set up directories
fprintf('setting up directories \n')
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
first_vol_filename = fullfile(data_root, 'first_vol_mamut');
data_filename = fullfile(data_root, 'mamut_data');
save_directory = fullfile(D, 'individual_neurons');
make_directory(save_directory);

% Initialize feature annotator
global feature_annotator_data;
feature_annotator_data.segmentation_options = struct();


%% Trim movie and append subsequent trials
fprintf('trimming movie and appending subsequent trials \n')
mov_index = data_file(end-5);
mov_name = @(x) sprintf('%s_0%i (%i).mat', strain, animal, 1+x);
mfile = matfile(data_file);
raw_buffer = mfile.images;
raw_buffer_trimmed = trim_movie(raw_buffer);
r0 = max_intensity_z(get_slice(raw_buffer_trimmed, 1));
num_movies = 5;
num_row = ceil(sqrt(num_movies-1));

% Append subsequent trials
i = 1;
while exist(mov_name(i), 'file') && i<num_movies
    data_file = fullfile(root, mov_name(i));
    mfile = matfile(data_file);
    raw_buffer = mfile.images;
    
    new_trimmed = trim_movie(raw_buffer);
    
    %align new_movie
    r1 = max_intensity_z(get_slice(new_trimmed, 1));
    figure(1); subplot(num_row, num_row*2, 2*i-1); imshowpair(r1, r0);

    [optimizer, metric] = imregconfig('multimodal');
    optimizer.InitialRadius = optimizer.InitialRadius/3.5; % Reduce the step size
    tform = imregtform(r1, r0, 'rigid', optimizer, metric);
    rotate_new = @(x) imwarp(x, tform, 'OutputView', imref2d(size(r0)));
    r1_transformed = rotate_new(r1);
    figure(1); subplot(num_row, num_row*2, 2*i); imshowpair(r1_transformed, r0);
    
    new_trimmed_rotated = rotate_new(new_trimmed);
    
    raw_buffer_trimmed = cat(4, raw_buffer_trimmed, ...
        new_trimmed_rotated);
    
    i = i+1;
%     pause(2);
end
num_trials = i;
raw_buffer_trimmed = LazyArray(raw_buffer_trimmed, @(x)x);
num_frames = size(raw_buffer_trimmed, 4);


%% Align red and green channels
fprintf('aligning red and green channels \n')

% Reference time for global rotation alignment.
t0 = 1;

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

r0 = get_slice(raw_red_z, t0);
g0 = get_slice(raw_green_z, t0);
figure(2); clf; subplot(121); imshowpair(g0, r0);

[optimizer, metric] = imregconfig('multimodal');
% optimizer.InitialRadius = optimizer.InitialRadius/3.5;
tform = imregtform(r0, g0, 'translation', optimizer, metric); tform_rg = tform;
rotate_red = @(x) imwarp(x, tform, 'OutputView', imref2d(size(r0)));

r0_transformed = rotate_red(r0);
figure(2); subplot(122); imshowpair(g0, r0_transformed);

registration_filename = fullfile(assets, 'red_green_alignment.png');
saveas(gcf, registration_filename);

green = raw_green;
red = LazyArray(raw_red, rotate_red);

green_scaled = RangeScaledArray(green);
red_scaled = RangeScaledArray(red);


%% Align to previous movie
fprintf('determining orientation \n')

data = ColorArray({red_scaled, green_scaled}, 4, eye(2), false);
data_rgb = ColorArray({red_scaled, green_scaled}, 4, eye(3), true);

%apply guess rotation from orienta_worm.m
transformation_location = fullfile(assets, 'transformation');
load(transformation_location); % load flipped and t_rot
data_rotated = ImwarpedArray(data);
data_rotated_rgb = ImwarpedArray(data_rgb);
data_rotated_test = ImwarpedArray.update_tform(data_rotated_rgb, 1, t_rot);
if flipped
    data_rotated_test = LazyArray(data_rotated_test, @flip);
end
data_rotated_test_z = LazyArray(data_rotated_test, @max_intensity_z);
figure(4); clf; imagesc(get_slice(data_rotated_test_z, 1));

%fine-tune alignment
target_stack = load_tiff_stack(D, 1);
target = max_intensity_z(target_stack); 
source = get_slice(get_slice(data_rotated_test_z,1),1);
[optimizer, metric] = imregconfig('multimodal');
tform = imregtform(source, target, 'rigid', optimizer, metric);
rotate_source = @(x) imwarp(x, tform, 'OutputView', imref2d(size(target)));
guess = rotate_source(source);
figure(2);  imshowpair(target, guess);


%% If it looks good, apply transformation to all frames
% otherwise, repeat previous cell
fprintf('orienting all frames \n')

for i = 1:num_frames
    data_rotated = ImwarpedArray.update_tform(data_rotated, i, t_rot);
    data_rotated_rgb = ImwarpedArray.update_tform(data_rotated_rgb, ...
        i, t_rot);
end

if flipped
    data_rotated = LazyArray(data_rotated, @flip);
    data_rotated_rgb = LazyArray(data_rotated_rgb, @flip);
end

view = imref2d(size(target));
data_aligned = ImwarpedArray(data_rotated, tform, view);
data_aligned_rgb = ImwarpedArray(data_rotated_rgb, tform, view);


%% Filter and normalize
fprintf('filtering and normalizing \n')

bump_finder = bump_detector([1,1], [3,3], [51 51]);
normalization = 0.7 * sum_all(bump_finder(bump_finder>0)*255);

dc = @(x) deconvlucy(x, fspecial('gaussian', 5, 2), 1);
flt = @(x, h) imfilter(double(dc(x)), h, 'circular');
flt_and_normalize = @(x) uint8(flt(x, bump_finder)/normalization*256-5);

% flt_normalize_slice = @(x) flt_and_normalize(x(crop_frames,...
%     :,:,:,:));

data_aligned_filtered_normalized = LazyArray(data_aligned, ...
    flt_and_normalize);
data_aligned_filtered_normalized_rgb = LazyArray(data_aligned_rgb, ...
    flt_and_normalize);

slice_num = [201 299];
figure(4); 
subplot(1,2,1); imshow(max_intensity_z(get_slice(data_aligned_rgb, slice_num(1))));
title(['Frame ' num2str(slice_num(2)) ' pre filt/deconv']);
subplot(1,2,2); imshow(max_intensity_z(get_slice(data_aligned_rgb, slice_num(2))));
title(['Frame ' num2str(slice_num(2)) ' pre filt/deconv']);
figure(5); 
subplot(1,2,1); imshow(max_intensity_z(get_slice(...
    data_aligned_filtered_normalized_rgb, slice_num(1))));
title(['Frame ' num2str(slice_num(1)) ' post filt/deconv']);
subplot(1,2,2); imshow(max_intensity_z(get_slice(...
    data_aligned_filtered_normalized_rgb, slice_num(2))));
title(['Frame ' num2str(slice_num(2)) ' post filt/deconv']);


%% Write out
% Estimated time ~1 hour. Run this, then begin proofreading annotation
fprintf('writing out \n')

tic;
% Write TIFs
red_stack = LazyArray(data_aligned_filtered_normalized,...
    @(x) squeeze(x(:,:,:,1)));
green_stack = LazyArray(data_aligned_filtered_normalized,...
    @(x) squeeze(x(:,:,:,2)));
array_to_TIFF_stacks(red_stack, red_stack_tiff_directory);
array_to_TIFF_stacks(green_stack, green_stack_tiff_directory);

% Save to big data viewer array
data_for_mamut = LazyArray(data_aligned_filtered_normalized,...
    @(x) permute(x(:,:,:,:,:), [2,1,3,4,5]));
write_big_data_viewer_array(data_filename, data_for_mamut);


%% Import tracks from MaMuT

tracks = tracks_from_mamut(data_filename, [1 1 2]);
features = features_from_tracks(tracks, num_frames);
feature_annotator_data.features = features;
save_features(features, D);


%% Prepare for tracking
v1 = load_tiff_stack(D, 1);
N_features = length(features);
cutoff = 0.80;
max_tries = 5;
discount = 0.95;
t_0 = 1;
track_size = [25 25 10];

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
fprintf('tracking kicks off \n');
tic;
parfor i = 1:length(features);
    warning off;
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
        modified_max_tries = max_tries;
        q = 0;
        f0 = f;
        while q<cutoff && ~isempty(q_remaining) && ...
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

        fprintf('Feature %d, time %d okay in %d tries. Score: %d. Parent: %d\n', ...
            i, t, tries, q, t_ref_best);
        
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
toc;
fprintf('tracking finished. \n');
feature_annotator_data.features = features;
save_features(features, D);



% %% Translate back to MaMuT
% 
% save_tracks = tracks_from_features(features);
% mamut_xml_from_tracks(save_tracks, data_filename, [1 1 2]);
% 
% 
% %% Extract calcium traces
% % [Optional] reads intensities of each channel at each time point, and
% % saves them as time series
% 
% feature_size = [3, 3, 4];
% options = struct(...
%     'npoints', 10 ...
% );
% 
% features = load_features(D);
% features = set_feature_size(features, feature_size);
% 
% rfp = get_intensity(features, red_stack_tiff_directory, options);
% gcamp = get_intensity(features, green_stack_tiff_directory, options);
% 
% 
% filename = fullfile(data_root, 'traces.mat');
% save(filename, 'gcamp', 'rfp');
% 
% figure(1); clf; plot(gcamp);
% saveas(gcf, fullfile(data_root, 'gcamp.jpg'));
% figure(2); clf; imagesc(transpose(gcamp));
% saveas(gcf, fullfile(data_root, 'heatmap.jpg'));
% 
% 
% 
% 
