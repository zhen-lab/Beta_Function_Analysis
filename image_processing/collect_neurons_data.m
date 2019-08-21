function [GCaMP, RFP, ratio, ratio_norm] = collect_neurons_data(folder, frames)

f = dir(fullfile(folder, '*.mat'));
numsamples = numel(f);
GCaMP = zeros(frames, numsamples);
RFP = zeros(frames, numsamples);
ratio = zeros(frames, numsamples);

for idx = 1 : numsamples
    name = fullfile(folder, f(idx).name);
    load(name, 'signal', 'signal_mirror');
    frame_length = min(size(signal{1,1},1), frames); % in case there are recordings with shorter time window
    GCaMP(1:frame_length, idx) = signal{1,1}(1:frame_length);
    RFP(1:frame_length, idx) = signal_mirror{1,1}(1:frame_length);
    ratio(1:frame_length, idx) = signal{1,1}(1:frame_length)./signal_mirror{1,1}(1:frame_length);
end

ratio_start = repmat(ratio(1, :), frames, 1);
ratio_norm = (ratio./ratio_start)';

GCaMP = GCaMP'; RFP = RFP'; ratio = ratio';

close all;
plot(ratio_norm');

end