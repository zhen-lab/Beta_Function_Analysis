function features = features_from_tracks(tracks, size_T)
% features = FEATURES_FROM_TRACKS(tracks)
%
%   Converts tracks (e.g. from TRACKS_FROM_MAMUT) to features for
%   processing.
%
% Author: Vivek Venkatachalam (vivekv2@gmail.com)

filter = @(f) tracks(f(tracks),:);

all_skeletonIds = unique(tracks.skeletonId);
%valid_SIDs = all_skeletonIds(all_skeletonIds>-1);
valid_SIDs = all_skeletonIds;

if length(valid_SIDs) > 1
    neuron_names = valid_SIDs;
    annotation = false;
else
    neuron_names = unique(tracks.id);
    annotation = true;
end


feature_size = [10 10 3];

features = {};
for i = 1:length(neuron_names)
    sid = neuron_names(i);
    if annotation
        track = filter(@(x) x.id == sid);
    else
        track = filter(@(x) x.skeletonId == sid);
    end
    
    f = new_feature(num2str(sid), size_T, feature_size);
    
    for j = 1:size(track,1)
        t = track{j, 'time'};
        x = track{j, 'x'};
        y = track{j, 'y'};
        z = track{j, 'z'};
        confidence = track{j, 'confidence'};
        f.coordinates(t,:) = [x y z];
        if confidence == -1
            f.is_registered(t) = true;
        end
    end
    
    f = translate_feature(f, -f.ref_offset);
    
    features{i} = f;
end