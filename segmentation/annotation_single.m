function [] = annotation_single(strain, animal)

%% Set up directories
tic;
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

%% Annotate neurons
Dirc = dir([D, '\*.tif']);
num_frames = length(Dirc(not([Dirc.isdir])));
v = load_tiff_stack(D, 1);
fprintf('interpolating volume data \n');
vq = interp3(double(v)); % interpolation of data to avoid consecutive maxima

fprintf('finding local maxima \n');
% num_neurons_to_identify = 250;
% maxima = find_local_maxima(v, [5, 5, 3], 'max', num_neurons_to_identify, ...
%     'minimum_separation', [0.5, 0.5, 0.5]);
maxima = find_local_maxima(vq, [5, 5, 3]);
maxima = ceil(maxima/2);
flength = size(maxima, 1);
fsize = [50, 50, 5];

features = cell(1, flength);
for i = 1:flength
    name = sprintf('neuron_%03d', i);
    f = new_feature(name, num_frames, fsize);
    f = translate_feature(f, maxima(i,:) - fsize/2);
    f.modified_coordinates(1,:) = f.coordinates(1,:);
    features{i} = f;
end

fprintf('saving data \n');
feature_annotator_data.features = features;
save_features(features, D);

save_tracks = tracks_from_features(features);
mamut_xml_from_tracks(save_tracks, data_filename, [1 1 2]);

toc;

end