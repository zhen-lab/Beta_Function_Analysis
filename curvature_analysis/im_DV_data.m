function [s] = im_DV_data(dorsal_ratio, ventral_ratio, range_overlap, range_separated, choice)

close all;

color_ctrl = [0.8 0.8 0.8; 0 0 0];
color_atr = [0.8 0.8 0.8; 0 0 0];
% color_atr = [0.5 1 0.5; 0 0.9 0];
clr = color_ctrl*choice(1) + color_atr*choice(2);

dorsal_ratio(dorsal_ratio==0) = NaN;
ventral_ratio(ventral_ratio==0) = NaN;

s(1) = figure(1);hold on
plot(dorsal_ratio', 'color', clr(1,:), 'linewidth', 2);
ylim(range_overlap);
plot(nanmean(dorsal_ratio), 'color', clr(2,:), 'linewidth', 4);
% plot([1 length(dorsal_ratio)], [1 1], '--', 'color', 0.9*[1 1 1], 'linewidth', 2);
% plot([1 length(dorsal_ratio)], [1 1], '--', 'color', [.5 1 .5], 'linewidth', 2);
% a = ylim;
% plot([0 0],[a(1) a(2)],'color','k');
set(gca, 'visible', 'off')
set(gcf, 'color', 'w')
% c = colorbar;
% set(c, 'visible', 'off');

s(2) = figure(2);hold on
plot(ventral_ratio', 'color', clr(1,:), 'linewidth', 2);
ylim(range_overlap)
plot(nanmean(ventral_ratio), 'color', clr(2,:), 'linewidth', 4);
% plot([1 length(dorsal_ratio)], [1 1], '--', 'color', 0.9*[1 1 1], 'linewidth', 2);
% plot([1 length(dorsal_ratio)], [1 1], '--', 'color', [.5 1 .5], 'linewidth', 2);
set(gca, 'visible', 'off')
set(gcf, 'color', 'w')
% c = colorbar;
% set(c, 'visible', 'off');

s(3) = figure(3);
imagesc(dorsal_ratio);
caxis(range_separated);
set(gca, 'ticklength', [0 0], 'xticklabel', [], 'yticklabel', []);
c = colorbar; colormap(flipud(gray));
set(c, 'yticklabel', []);

s(4) = figure(4);
imagesc(ventral_ratio);
caxis(range_separated);
set(gca, 'ticklength', [0 0], 'xticklabel', [], 'yticklabel', []);
c = colorbar; colormap(flipud(gray));
set(c, 'yticklabel', []);

end