function [s] = im_DV_data_overlaid(input_d, input_v, range_overlap, lw, wd)

close all;

dnum = size(input_d,1);
vnum = size(input_v,1);

if dnum~=2 || vnum~=2
    
    fprintf('dimension of input files are incorrect.\n');

else
    
    color_ctrl = [0.8 0.8 0.8; 0 0 0];
    color_exp = [0.95 0.6 0.95; 0.95 0 0.95];

    input_d_exp = input_d{1}; input_d_exp(input_d_exp==0) = NaN;
    input_d_ctrl = input_d{2}; input_d_ctrl(input_d_ctrl==0) = NaN;
    input_v_exp = input_v{1}; input_v_exp(input_v_exp==0) = NaN;
    input_v_ctrl = input_v{2}; input_v_ctrl(input_v_ctrl==0) = NaN;
    
    s(1) = figure(1); % Dorsal
    hold on;
    plot(input_d_ctrl', 'color', color_ctrl(1,:), 'linewidth', lw);
    plot(nanmean(input_d_ctrl), 'color', color_ctrl(2,:), 'linewidth', 3*lw);
    plot(input_d_exp', 'color', color_exp(1,:), 'linewidth', lw);
    plot(nanmean(input_d_exp), 'color', color_exp(2,:), 'linewidth', 3*lw);
    ylim(range_overlap);
    set(gca, 'visible', 'off');
    set(gcf, 'color', 'w', 'units', 'normalized', 'outerposition', [0 0 wd 1]);

    s(2) = figure(2); % Ventral
    hold on;
    plot(input_v_ctrl', 'color', color_ctrl(1,:), 'linewidth', lw);
    plot(nanmean(input_v_ctrl), 'color', color_ctrl(2,:), 'linewidth', 3*lw);
    plot(input_v_exp', 'color', color_exp(1,:), 'linewidth', lw);
    plot(nanmean(input_v_exp), 'color', color_exp(2,:), 'linewidth', 3*lw);
    ylim(range_overlap);
    set(gca, 'visible', 'off');
    set(gcf, 'color', 'w', 'units', 'normalized', 'outerposition', [0 0 wd 1]);
    
end

end