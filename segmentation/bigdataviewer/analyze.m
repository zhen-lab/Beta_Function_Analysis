%%
root = 'V:\data\170619 BigDataViewer compatibility';
cd(root);
addpath(genpath(fullfile(root, 'src')));

dataset_name = 'ZM9624_OH11119_FYNJC_109-v1';
mamut_file = fullfile(root, dataset_name);

tracks = tracks_from_mamut(mamut_file);

%%
filter_data = @(d,f) d(f(d),:);

sample_track_1 = filter_data(tracks, @(x) x.skeletonId == 60);
sample_track_2 = filter_data(tracks, @(x) x.skeletonId == 82);
sample_track_3 = filter_data(tracks, @(x) x.skeletonId == 112);

sample_tracks = [sample_track_1; sample_track_2; sample_track_3];

features = features_from_tracks(sample_tracks);