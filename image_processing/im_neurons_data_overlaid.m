function [s] = im_neurons_data_overlaid(ratio_norm, range_overlap, lw)

close all;

numgroup = size(ratio_norm,2);
switch numgroup
    case 1
        color_exp = [0.95 0.6 0.95; 0.95 0 0.95];
        ratio_norm_exp = ratio_norm{:}; ratio_norm_exp(ratio_norm_exp==0) = NaN;
        s = figure; hold on;
        plot(ratio_norm_exp', 'color', color_exp(1,:), 'linewidth', lw);
        plot(nanmean(ratio_norm_exp), 'color', color_exp(2,:), 'linewidth', 3*lw);
        ylim(range_overlap);
        set(gca, 'visible', 'off')
        set(gcf, 'color', 'w')
    case 2
        color_ctrl = [0.8 0.8 0.8; 0 0 0];
        color_exp = [0.95 0.6 0.95; 0.95 0 0.95];
        % color_atr = [0.5 1 0.5; 0 0.9 0];
        % clr = color_ctrl*choice(1) + color_atr*choice(2);
        ratio_norm_exp = ratio_norm{1}; ratio_norm_exp(ratio_norm_exp==0) = NaN;
        ratio_norm_ctrl = ratio_norm{2}; ratio_norm_ctrl(ratio_norm_ctrl==0) = NaN;
        s = figure; hold on;
        plot(ratio_norm_ctrl', 'color', color_ctrl(1,:), 'linewidth', lw);
        plot(nanmean(ratio_norm_ctrl), 'color', color_ctrl(2,:), 'linewidth', 3*lw);
        plot(ratio_norm_exp', 'color', color_exp(1,:), 'linewidth', lw);
        plot(nanmean(ratio_norm_exp), 'color', color_exp(2,:), 'linewidth', 3*lw);
        ylim(range_overlap);
        set(gca, 'visible', 'off')
        set(gcf, 'color', 'w')
end

end