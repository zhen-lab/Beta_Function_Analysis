function propagate_tracking_traces(strain, animal)

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

% Load features.mat
tic;
load([D '\features.mat']);
toc;
%% Translate back to MaMuT
tic;
fprintf('saving mamut file \n');
save_tracks = tracks_from_features(features);
mamut_xml_from_tracks(save_tracks, data_filename, [1 1 2]);
toc;

%% Extract calcium traces
% [Optional] reads intensities of each channel at each time point, and
% saves them as time series
tic;
fprintf('extracting calcium traces \n');
feature_size = [3, 3, 4];
options = struct(...
    'npoints', 10 ...
);

features = load_features(D);
features = set_feature_size(features, feature_size);

rfp = get_intensity(features, red_stack_tiff_directory, options);
gcamp = get_intensity(features, green_stack_tiff_directory, options);

fprintf('saving calcium traces \n');
filename = fullfile(data_root, 'traces.mat');
save(filename, 'gcamp', 'rfp');

figure(1); clf; plot(gcamp);
saveas(gcf, fullfile(data_root, 'gcamp.jpg'));
figure(2); clf; imagesc(transpose(gcamp));
saveas(gcf, fullfile(data_root, 'heatmap.jpg'));
toc;
end
