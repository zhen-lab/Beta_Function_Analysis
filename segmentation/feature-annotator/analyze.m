%%


root = 'D:\ZM9624_temp_timelapse';
file_pattern = '1213_ZM9624_001.mat';

interesting_frames = [171:183, 232:244, 524:534, 994:1014];
valid_frames = [1:1200];

run_number = 0928;
data_file = fullfile(root,file_pattern);

mfile = matfile(data_file);

% Set a location for processed data.
[~, f, ~] = fileparts(data_file);
data_root = fullfile(root, f);
make_directory(data_root);

% Create a subfolder to hold useful summary data.
assets = fullfile(data_root, 'assets');
make_directory(assets);

% Create folders for image stacks
red_filtered = fullfile(data_root, 'red_filtered');
green_filtered = fullfile(data_root, 'green_filtered');
make_directory(red_filtered);
make_directory(green_filtered);

cd(root);

feature_annotator();

global feature_annotator_data;
feature_annotator_data.segmentation_options = struct();

%%

raw_buffer = MatfileArray(data_file, 'images');
raw_buffer_cropped = TimeMappedArray(raw_buffer, valid_frames);

%% Create an aligned array

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

raw_image = LazyArray(raw_buffer_cropped, get_zyla_vol);
raw_green = LazyArray(raw_image, get_green_slices);
raw_red = LazyArray(raw_image, get_red_slices);

raw_green_z = LazyArray(raw_green, @max_intensity_z);
raw_red_z = LazyArray(raw_red, @max_intensity_z);

r0 = get_slice(raw_red_z, t0);
g0 = get_slice(raw_green_z, t0);
figure(2); clf; subplot(121); imshowpair(g0, r0);

[optimizer, metric] = imregconfig('monomodal');
tform = imregtform(r0, g0, 'rigid', optimizer, metric);
rotate_red = @(x) imwarp(x, tform, 'OutputView', imref2d(size(r0)));

r0_transformed = rotate_red(r0);
figure(2); subplot(122); imshowpair(g0, r0_transformed);

registration_filename = fullfile(assets, 'red_green_alignment.png');
saveas(gcf, registration_filename);

green = raw_green;
red = LazyArray(raw_red, rotate_red);

%% @Jasper: make two arrays called 'red' and 'green', each 4D.
% Hopefully the rest of this will work if you do that...


%% Prepare a combined array with the correct scale along with MIPs

green_scaled = RangeScaledArray(green);
red_scaled = RangeScaledArray(red);

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

mips = ColorArray({red_mips, green_mips}, 3, [1,0; 0,1; 0,0], true);
mips = AnnotatedArray(mips);
mips.Functions{1} = @AnnotatedArray.timestamp;
mips.FunctionArgs{1} = {[10, 10]};

figure(3); clf; imshow(get_slice(mips_side_by_side, 10));

mips = ColorArray({red_mips, green_mips}, 3, eye(3), true);
figure(4); clf; imshow(autoscale(get_rgb_slice(mips, 1)));
figure(5); clf; imshow(autoscale(get_rgb_slice(mips, size(mips, 4))));

%% Write the MIPs to a video.

mipsfile_1 = fullfile(assets, 'MIP_side_by_side.mp4');
video_from_array(mips_side_by_side, mipsfile_1, ...
    'framerate', 10);

mipsfile_2 = fullfile(assets, 'MIP.mp4');
video_from_array(mips, mipsfile_2, ...
    'framerate', 10);

%% Determine a canonical orientation by clustering.
tic

smoothing = 4;
decimation = 2;

N_sample = 100;
N_tot = size(red_scaled, 4);

red_mip_z = LazyArray(red_scaled, @max_intensity_z);

smth = @(x,y) automatically_smooth_image(x, y);
shrink = @(x,y) x(1:y:end, 1:y:end, :);
dec = @(x, y, z) shrink(smth(x,y),z);

thumbs = LazyArray(red_mip_z, @(x) dec(x, smoothing, decimation));

A = CachedArray(LazyArray(thumbs, @center_image));

sample_times = 1:floor(N_tot/N_sample):N_tot;

N = length(sample_times);
C = zeros(N, N);

for i = 1:N
    
    for j = (i+1):N
        C(i,j) = get_image_overlap(get_slice(A,i), get_slice(A,j));
    end
end

C = C + C' + eye(N);
toc

%%

Y = 1-C;
figure(1); clf; imagesc(Y);

score = sum(Y);
[~, order] = sort(score);

target_idx = order(1);
target_time = sample_times(target_idx);

target = get_slice(A, target_time);

figure(2); imshow(target);

%%

for i = 1:100
    figure(2); clf; imshow(get_slice(A, order(i))*5);
    drawnow;
end

%% update target angle

target = get_slice(A, target_time);
additional_angle = 30;
target = imrotate(target, additional_angle, 'nearest', 'crop');
figure(2); imagesc(target);

%% Determine rigid transformation

get_im = @(x) get_slice(thumbs, x);

siz = size(thumbs);

center = siz(1:2)/2;
t_init = affine2d(eye(3));
initial_angle = 180;

t_init = rotate_tform(t_init, initial_angle, center);
t_last = t_init;

rigid_thumbs = ImwarpedArray(thumbs);

fits = [];
for i = valid_frames
    
    im = get_slice(thumbs,i);
    
    [new, t_last] = register_images(im, target, 'guess', t_last, ...
        'correlation_target', 0.75);
    
    t_last = rotate_tform(t_last, additional_angle, center);
    
    rigid_thumbs = ImwarpedArray.update_tform(rigid_thumbs, i, t_last);
    
    fit = get_image_overlap(new, target);
    fits(i) = fit;
    
    title([num2str(i) ' fit: ' num2str(fit)]);
    pause(0.01);
    
end

 %% Save thumbnails.

save(fullfile(data_root, 'rigid_thumbs.mat'), 'rigid_thumbs');

%% Convert rigid transformations to larger view.

data_rgb = ColorArray({red_scaled, green_scaled}, 4, eye(3), true);
data = ColorArray({red_scaled, green_scaled}, 4, eye(2), false);

data_aligned = ImwarpedArray.move_tforms(rigid_thumbs, data);
data_aligned_rgb = ImwarpedArray.move_tforms(rigid_thumbs, data_rgb);

data_aligned_rgb_z = LazyArray(data_aligned_rgb, @max_intensity_z);

for i = 1:10:1200
    im = get_slice(data_aligned_rgb_z, i);
    figure(3); imagesc(im);
    drawnow;
end

%% Trim edges

data_xy = LazyArray(data_aligned, @(x) max_intensity_z(...
    max_intensity(x,4)));

S = size(data_xy);
N_frames = S(end);
sample_frames = 1:floor(N_frames/50):N_frames;
data_xy_sampled = TimeMappedArray(data_xy, sample_frames);

[~, idx] = trim_edges(get_array_data(data_xy_sampled), ...
    'filter', @(x) imfilter(x-2, ones(4,4)/16));

data_trimmed = LazyArray(data_aligned, @(x) x(idx{1}, idx{2}, :, :, :));
data_trimmed_rgb = LazyArray(data_aligned_rgb, ...
    @(x) x(idx{1}, idx{2}, :, :, :));

data_trimmed_rgb_z = LazyArray(data_trimmed_rgb, @max_intensity_z);

%% Check trimming

figure(1); image(max_intensity_z(get_slice(data_trimmed_rgb, 200)));

for i = 1:size(data_trimmed_rgb_z, 4)
    imagesc(get_slice(data_trimmed_rgb_z, i));
    pause(0.01);
    drawnow;
end

%% Make movies

mips = ThreeViewArray(data_trimmed_rgb, pixel_ratio);

mips = AnnotatedArray(mips);
mips.Functions{1} = @AnnotatedArray.timestamp;
mips.FunctionArgs{1} = {[10, 10]};

mipsfile_3 = fullfile(assets, 'MIPs_aligned.mp4');
video_from_array(mips, mipsfile_3, ...
    'framerate', 10);

%% Save array (link)

filename = fullfile(data_root, 'data_aligned_and_trimmed.mat');
save(filename, 'data_trimmed');


%% Load array

filename = fullfile(data_root, 'data_aligned_and_trimmed.mat');
s = load(filename);
data_trimmed = s.data_trimmed;

%% Add a filter

size_data_trimmed = size(data_trimmed);

idx{1} = 1:256;
idx{2} = 1:174;
idx{3} = 1:20;

h = fspecial('gaussian', 5, 1);

h_full = @(x) imfilter(x(idx{1}, idx{2}, idx{3}, :), h);

data_trimmed_and_filtered = LazyArray(data_trimmed, h_full);


%% Write out

red_stack_tiff_directory = fullfile(data_root, 'red');
green_stack_tiff_directory = fullfile(data_root, 'green');
make_directory(red_stack_tiff_directory);
make_directory(green_stack_tiff_directory);

red_stack = LazyArray(data_trimmed, @(x) squeeze(x(:,:,:,1)));
green_stack = LazyArray(data_trimmed, @(x) squeeze(x(:,:,:,2)));

array_to_TIFF_stacks(red_stack, red_stack_tiff_directory);
array_to_TIFF_stacks(green_stack, green_stack_tiff_directory);

%% Write out filtered

red_stack = LazyArray(data_trimmed_and_filtered, ...
    @(x) squeeze(x(:,:,:,1)));
green_stack = LazyArray(data_trimmed_and_filtered, ...
    @(x) squeeze(x(:,:,:,2)));

array_to_TIFF_stacks(red_stack, red_filtered);
array_to_TIFF_stacks(green_stack, green_filtered);

%%
D = red_stack_tiff_directory;
v1 = load_tiff_stack(D, 1);

%%

maxima = find_local_maxima(v1, [5, 5, 3], 'max', 50, ...
    'minimum_separation', [3, 3, 1]);

size_T = get_size_T(D);
fsize = [50, 50, 5];

features = {};
for i = 1:size(maxima, 1)
    
    name = sprintf('neuron_%03d', i);
    f = new_feature(name, size_T, fsize);
    f = translate_feature(f, maxima(i,:) - fsize/2);
    f.modified_coordinates(1,:) = f.coordinates(1,:);
    
    features{i} = f;
    
end

feature_annotator_data.features = features;
save_features(features, red_filtered);

%% De-duplicate and rename

separation = [3,3,1];

name_fmt = 'neuron_%03d';

features = load_features(red_filtered);

i = 1;
while i < length(features)
    
    c0 = features{i}.modified_coordinates(1,:);
    
    marked = [];
    for j = (i+1):length(features)
        
        cj = features{j}.modified_coordinates(1,:);
        if all(abs(cj-c0) <= separation)
            marked = [marked j];
        end
        
    end
    features(marked) = [];
    
    i = i + 1;
end


for i = 1:length(features) 
    features{i}.name = sprintf(name_fmt, i);
end

save_features(features, red_filtered);

%% Get rid of registration information.

features = load_features(red_stack_tiff_directory);

for i = 1:length(features)
    
    features{i}.registration_info = [];
    
end

save_features(features, red_stack_tiff_directory);

%% 

D = red_filtered;
save_directory = fullfile(D, 'individual_neurons');
make_directory(save_directory);
cutoff = 0.1;
max_tries = 30;
discount = 0.95;

track_size = [35 35 6];

features = load_features(D);

features = set_feature_size(features, track_size);

feature_annotator_data.features = features;
feature_annotator_data.image_location = D;

N_features = length(features);
size_T = get_size_T(D);

t_0 = 1;

options = struct(...
    'drift', [2, 2, 1], ...
    'imreg_padding', [10, 10, 2], ...
    't_ref', t_0 ...
);

feature_annotator_data.segmentation_options = ...
    merge_struct(options, feature_annotator_data.segmentation_options);

%% @Jasper: The stuff below here is a bit stale.
% Go back to using feature_annotator (or writing your own registration
% loops) at this point.


%%
parfor i = 1:N_features;

    f = features{i};

    disp(['Registering feature ', num2str(i)]);

    good_frames = get_modified_frames(f);

    qualities = zeros(1, size_T);
    qualities(good_frames) = 1;
    
    depths = ones(1, size_T)*1000;
    depths(good_frames) = 0;

    for t = 1:size_T-1%size_T-1

        t_ref = t_0;
        
        q_remaining = qualities .* discount.^(depths);
        t_remaining = 1:size_T;

        tries = 0;
        modified_max_tries = min([t-1, max_tries]);
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

        disp(sprintf('Time %d okay in %d tries. Score: %d. Parent: %d', ...
            t, tries, q, t_ref_best));
        
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

save_features(features, D);

%%

D = red_filtered;
save_directory = fullfile(D, 'individual_neurons');
make_directory(save_directory);
cutoff = 0.75;
max_tries = 30;
discount = 0.95;

track_size = [35 35 4];

features = load_features(D);

features = set_feature_size(features, track_size);
feature_annotator_data.features = features;

N_features = length(features);
size_T = get_size_T(D);

t_0 = 1;

options = struct(...
    'drift', [2, 2, 1], ...
    'imreg_padding', [5, 5, 1], ...
    't_ref', t_0 ...
);

feature_annotator_data.segmentation_options = merge_struct(options, ...
    feature_annotator_data.segmentation_options);

%%
warning('off', 'images:regmex:registrationOutBoundsTermination');

%%
parfor i = 17:N_features;

    f = features{i};

    disp(['Registering feature ', num2str(i)]);

    good_frames = get_modified_frames(f);

    qualities = zeros(1, size_T);
    qualities(good_frames) = 1;
    
    depths = ones(1, size_T)*1000;
    depths(good_frames) = 0;

    for t = 2:size_T-1

        t_ref = t_0;
        
        q_remaining = qualities .* discount.^(depths);
        t_remaining = 1:size_T;

        tries = 0;
        modified_max_tries = min([t-1, max_tries]);
        q = 0;
        f0 = f;
        while q<cutoff && ~isempty(q_remaining) ...
           && tries < modified_max_tries

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

        disp(sprintf('Time %d okay in %d tries. Score: %d. Parent: %d', ...
            t, tries, q, t_ref_best));
        
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

save_features(features, D);

%% Get features from individual files

features = load_features(red_filtered);
save_directory = fullfile(red_filtered, 'individual_neurons');

for i = 1:N_features
    filename = fullfile(save_directory, sprintf('neuron_%04d.mat', i));
    if exist(filename, 'file')
        n = load_features(filename);
        features{i} = n{1};
    end
end
    



%%

feature_size = [3, 3, 2];
options = struct(...
    'npoints', 10 ...
);

features = load_features(red_filtered);
features = set_feature_size(features, feature_size);

gcamp = get_intensity(features, green_stack_tiff_directory, options);
rfp = get_intensity(features, red_stack_tiff_directory, options);

filename = fullfile(data_root, 'traces.mat');
save(filename, 'gcamp', 'rfp');

%%

S = load(fullfile(data_root, 'traces.mat'));
rfp = S.rfp;
gcamp = S.gcamp;

low = 5;
hi = 100;

for i = 1:size(gcamp, 2)
    
    g = gcamp(:,i);
    
    g0 = causal_smooth_ignore_nans(g, hi);
    g1 = g - g0;
    g2 = causal_smooth_ignore_nans(g1, low);
    
    s(:,i) = g2;
    
end

%%



%% Check registration

feature_idx = 59;
t_idx = 1:100;

feature = feature_annotator_data.features{feature_idx};

size_T = size(feature.coordinates, 1);

scores = zeros(size_T, 1);

threshold = 0.8;
good_fig = 1;
bad_fig = 2;

for t = t_idx
    
    s = feature.registration_info{t}.image_registration_score;
    scores(t) = s;
    
end

s = scores(t_idx);
x = feature.coordinates(t_idx, 2);
y = feature.coordinates(t_idx, 1);
z = feature.coordinates(t_idx, 3);

figure(1); clf; plot(scores, 'b.');
hold on;

%bad_idx = [37 38 42 44 47 97];
bad_idx = [27 33 35 41 44 45 47 84 86];
plot(bad_idx, x(bad_idx), 'r.', 'MarkerSize', 10);


figure(2); clf;
subplot(131);
plot(x, 'b.');
hold on;
plot(bad_idx, x(bad_idx), 'r.', 'MarkerSize', 10);
subplot(132);
plot(y, 'b.');
hold on;
plot(bad_idx, y(bad_idx), 'r.', 'MarkerSize', 10);
subplot(133);
plot(z, 'b.');
hold on;
plot(bad_idx, z(bad_idx), 'r.', 'MarkerSize', 10);

%% medfilt1 on coordinates to determine bad frames

