%% Prepare variables if workspace from previous sections is cleared
% Set up directories

% USER: Set these to directory of dataset
root = pwd; %folder containing dataset 
file_pattern = 'AVA_01 (3).mat'; %filename of dataset (e.g. AVA_01 (1).mat)

data_file = fullfile(root, file_pattern);
mfile = matfile(data_file);

[~, f, ~] = fileparts(data_file);
data_root = fullfile(root, f);
cd(root);

red_stack_tiff_directory = fullfile(data_root, 'red');
green_stack_tiff_directory = fullfile(data_root, 'green');
D = red_stack_tiff_directory;
data_filename = fullfile(data_root, 'mamut_data');
save_directory = fullfile(D, 'individual_neurons');
make_directory(save_directory);

feature_annotator_data.segmentation_options = struct();

num_frames = 200;

%% Import tracks from MaMuT

tracks = tracks_from_mamut(data_filename, [1 1 2]);
features = features_from_tracks(tracks, num_frames);

feature_annotator_data.features = features;
save_features(features, D);

disp('Features obtained');

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

disp('Tracking prepared');

%% Track neurons

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
        f0 = f;
        while q<cutoff && ~isempty(q_remaining) && ...
              tries < modified_max_tries

            idx_next = find(q_remaining == max(q_remaining), 1);
            disp(num2str(idx_next));
            
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

disp('Tracking finished');

%% Translate back to MaMuT

save_tracks = tracks_from_features(features);
mamut_xml_from_tracks(save_tracks, data_filename, [1 1 2]);

disp('Tracks saved');

%% NOW PROOFREAD IN MaMuT
%% Re-import proofread tracks from MaMuT

tracks = tracks_from_mamut(data_filename, [1 1 2]);
features = features_from_tracks(tracks, num_frames);

feature_annotator_data.features = features;
save_features(features, D);


%% Extract calcium traces
% [Optional] reads intensities of each channel at each time point, and
% saves them as time series

feature_size = [3, 3, 4];
options = struct(...
    'npoints', 10 ...
);

features = load_features(D);
features = set_feature_size(features, feature_size);

rfp = get_intensity(features, red_stack_tiff_directory, options);
gcamp = get_intensity(features, green_stack_tiff_directory, options);


filename = fullfile(data_root, 'traces.mat');
save(filename, 'gcamp', 'rfp');

figure(1); clf; plot(gcamp);
saveas(gcf, fullfile(data_root, 'gcamp.jpg'));
figure(2); clf; imagesc(transpose(gcamp));
saveas(gcf, fullfile(data_root, 'heatmap.jpg'));


%% Load features from files

features = load_features(D);

for i = 1:length(features)
    filename = fullfile(save_directory, sprintf('neuron_%04d.mat', i));
    if exist(filename, 'file')
        n = load_features(filename);
        features{i} = n{1};
    end
end
    