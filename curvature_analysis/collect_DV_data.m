function [dorsal_ratio, ventral_ratio, dorsal, ventral] = collect_DV_data(folder, frames)

f = dir(fullfile(folder, '*.mat'));
numsamples = numel(f);
dorsal = zeros(frames, numsamples);
ventral = zeros(frames, numsamples);
headneck = 36;

for idx = 1 : numsamples
    name = fullfile(folder, f(idx).name);
    load(name, 'dorsal_smd', 'dorsal_smd_r', 'ventral_smd', 'ventral_smd_r');
    frame_length = min(size(dorsal_smd,2), frames); % in case there are recordings with shorter time window
    dorsal(1:frame_length, idx) = mean(dorsal_smd(headneck:end, 1:frame_length)./dorsal_smd_r(headneck:end, 1:frame_length))';
    ventral(1:frame_length, idx) = mean(ventral_smd(headneck:end, 1:frame_length)./ventral_smd_r(headneck:end, 1:frame_length))';
end

dorsal_start = repmat(dorsal(1, :), frames, 1);
ventral_start = repmat(ventral(1, :), frames, 1);
dorsal_ratio = (dorsal./dorsal_start)';
ventral_ratio = (ventral./ventral_start)';

dorsal = dorsal'; ventral = ventral';

end