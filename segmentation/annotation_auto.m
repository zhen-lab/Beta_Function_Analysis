tic;

lst = dir; 
root = pwd;
num_neurons_to_identify = 130;

for j = 1:length(lst)-2
    
    file_pattern = lst(j+2).name; 
    D = fullfile(root, file_pattern, 'red');
    data_filename = fullfile(root, file_pattern, 'mamut_data');
    
    Dirc = dir([D, '\*.tif']);
    num_frames = length(Dirc(not([Dirc.isdir])));
    v = load_tiff_stack(D, 1);

    maxima = find_local_maxima(v, [5, 5, 3], 'max', ...
        num_neurons_to_identify, 'minimum_separation', [3, 3, 1]);
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

    feature_annotator_data.features = features;
    save_features(features, D);

    save_tracks = tracks_from_features(features);
    mamut_xml_from_tracks(save_tracks, data_filename, [1 1 2]);
    
end

toc;

clear;