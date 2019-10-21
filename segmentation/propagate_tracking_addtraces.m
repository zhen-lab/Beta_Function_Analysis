function propagate_tracking_addtraces(strain, animal)

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

tracks = tracks_from_mamut(data_filename, [1 1 2]); num_frames = 500;
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

end