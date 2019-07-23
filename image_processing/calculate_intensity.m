function [F_denoised, F_sorted, F_thresh] = calculate_intensity(F)

F_sorted = sort(F, 'descend');
thresh_signal = ceil(length(F) / 8);
thresh_baseline = floor(length(F) / 8);

F_signal = mean(F_sorted(1:thresh_signal)); % The highest 1/4 pixels are considered signals
F_baseline = mean(F_sorted((end - thresh_baseline):end)); % The lowest 1/8 pixels are considered baseline

F_thresh = F_sorted(thresh_signal);
F_denoised = F_signal - F_baseline; % average over the first n brightest pixels

end