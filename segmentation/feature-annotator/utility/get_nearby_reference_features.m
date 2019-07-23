function [nearest, nearest_features] = get_nearby_reference_features(...
    target_feature, all_features, t_ref, N, varargin)
% function [nearest, nearest_features] = get_nearby_reference_features(...
%   target_feature, all_reference_features, t_ref, N)
%
%   finds reference features near a given reference feature and returns
%   both the indices and the features themselves. t_ref is the time index
%   at which distances are calculated, and N is the number of features to
%   return;


default_options = struct(...
    'confidence', 0.95...
);
                
input_options = varargin2struct(varargin{:}); 
options = mergestruct(default_options, input_options);

reference_coords = [];
index = [];
all_reference_features = {};

confidence = options.confidence;

for n = 1:length(all_features)
    %make sure features are good
    conf = all_features{n}.is_registered(t_ref) > confidence...
        && all_features{n}.is_registered(t_ref - 1) > confidence;
    if n ~= target_feature && conf
        all_reference_features{end+1} = all_features{n};
        index(end+1) = n;
    end
end


for j = 1:length(all_reference_features)

    reference_coords(:,j) = ...
        get_feature_center(all_reference_features{j}, t_ref)';

end

feature_coords = get_feature_center(target_feature, t_ref)';
if length(reference_coords) > N
    nearest_indeces = find_nearest(feature_coords, reference_coords, N, 'scales', [1 1 5]);
    nearest = index(nearest_indeces);
    nearest_features = all_reference_features(nearest_indeces);
else
    disp('sorry, not enough confident features');
    return
end
