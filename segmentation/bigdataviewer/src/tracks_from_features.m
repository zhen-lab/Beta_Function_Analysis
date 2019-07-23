function tracks = tracks_from_features(features)

num_features = length(features);

if num_features > 0
    
    size_T = length(features{1}.coordinates);
    length_tracks = size_T*num_features;
    id =  transpose(1:length_tracks);
    type = zeros(length_tracks, 1)-1;
    y = NaN(length_tracks, 1);
    x = NaN(length_tracks, 1);
    z = NaN(length_tracks, 1);
    radius = 3*ones(length_tracks, 1);
    parent_id = zeros(length_tracks, 1);
    time = zeros(length_tracks, 1);
    confidence = 100*ones(length_tracks, 1);
    skeletonId = zeros(length_tracks, 1);
    
    for i = 1:num_features
        f = translate_feature(features{i}, features{i}.ref_offset);
        for j = 1:size_T
            index = (j-1)*num_features+i;
            x(index) = f.coordinates(j,1);
            y(index) = f.coordinates(j,2);
            z(index) = f.coordinates(j,3);
            time(index) = j;
            skeletonId(index) = i-1;
            if j > 1
                parent_id(index) = id(index-num_features);
            else
                parent_id(index) = -1;
            end
        end
    end
    
    tracks = table(id, type, y, x, z, radius, parent_id, ...
        time, confidence, skeletonId);
    
end
