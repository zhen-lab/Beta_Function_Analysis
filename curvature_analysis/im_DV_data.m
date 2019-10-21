function [] = im_DV_data(dorsal_ratio, ventral_ratio, range_overlap, range_separated, choice)

close all;

color_ctrl = [0.8 0.8 0.8; 0 0 0];
color_atr = [0.5 1 0.5; 0 0.9 0];
clr = color_ctrl*choice(1) + color_atr*choice(2);

figure(1);
plot(dorsal_ratio', 'color', clr(1,:), 'linewidth', 2);
ylim(range_overlap);
hold on
plot(mean(dorsal_ratio), 'color', clr(2,:), 'linewidth', 4);
set(gca, 'visible', 'off')
set(gcf, 'color', 'w')

figure(2);
plot(ventral_ratio', 'color', clr(1,:), 'linewidth', 2);
ylim(range_overlap)
hold on
plot(mean(ventral_ratio), 'color', clr(2,:), 'linewidth', 4);
set(gca, 'visible', 'off')
set(gcf, 'color', 'w')

figure(3);
imagesc(dorsal_ratio); colormap jet; 
caxis(range_separated);
set(gca, 'ticklength', [0 0], 'xticklabel', [], 'yticklabel', []);
c = colorbar;
set(c, 'yticklabel', []);

figure(4);
imagesc(ventral_ratio); colormap jet; 
caxis(range_separated);
set(gca, 'ticklength', [0 0], 'xticklabel', [], 'yticklabel', []);
c = colorbar;
set(c, 'yticklabel', []);

end