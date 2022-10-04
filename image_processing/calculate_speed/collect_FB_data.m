function [forward, backward, forward_pooled, backward_pooled] = collect_FB_data(folder, numneurons)

f = dir(fullfile(folder, '*.mat'));
numsamples = numel(f);
forward = {};
backward = {};
forward_pooled = cell(numsamples,numneurons);
backward_pooled = cell(numsamples,numneurons);

for idx = 1 : numsamples
    name = fullfile(folder, f(idx).name);
    load(name, 'signal_directions_separated_forward', 'signal_directions_separated_backward');
    for j = 1 : numneurons
        if ~isempty(signal_directions_separated_forward)
            forward_pooled{idx,j} = cat(1, signal_directions_separated_forward{:,j});
            forward = [forward; signal_directions_separated_forward];
        else
            forward_pooled{idx,j} = [];
            forward = [forward; cell(1,numneurons)];
        end
        if ~isempty(signal_directions_separated_backward)
            backward_pooled{idx,j} = cat(1, signal_directions_separated_backward{:,j});
            backward = [backward; signal_directions_separated_backward];        
        else
            backward_pooled{idx,j} = [];
            backward = [backward; cell(1,numneurons)];
        end
    end
end

forward_pooled_avg = cellfun(@nanmean, forward_pooled);
forward_pooled_avg_cat = cat(1, forward_pooled_avg(:));
backward_pooled_avg = cellfun(@nanmean, backward_pooled);
backward_pooled_avg_cat = cat(1, backward_pooled_avg(:));
pooled_avg_cat = [forward_pooled_avg_cat backward_pooled_avg_cat];

save('Forward_Backward_Pooled.mat', 'forward_pooled',  'backward_pooled', 'forward', 'backward', 'pooled_avg_cat');

end