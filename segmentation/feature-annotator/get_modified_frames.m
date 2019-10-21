function frames = get_modified_frames(feature)

mc = feature.modified_coordinates;

frames = [];
for i = 1:size(mc, 1)
    if all(~isnan(mc(i, :)))
        frames = [frames i];
    end
end