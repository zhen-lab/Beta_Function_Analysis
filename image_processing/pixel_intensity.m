function [F_denoised, F_sorted, F_thresh] = pixel_intensity(F, sig, base)

F_vec = reshape(F, [], 1);
F_sorted = sort(F_vec, 'descend');
thresh_signal = ceil(length(F) * sig);
thresh_baseline = floor(length(F) * base);

F_signal = mean(F_sorted(1:thresh_signal)); % The highest pixels are considered signals
F_baseline = mean(F_sorted((end - thresh_baseline):end)); % The lowest pixels are considered baseline

F_thresh = F_sorted(thresh_signal);
F_denoised = F_signal - F_baseline; % average over the first n brightest pixels

end