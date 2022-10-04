function [s] = im_neurons_data(ratio_norm, range_overlap, range_separated, choice)

close all;

color_ctrl = [0.8 0.8 0.8; 0 0 0];
color_atr = [0.8 0.8 0.8; 0 0 0];
% color_atr = [0.5 1 0.5; 0 0.9 0];
clr = color_ctrl*choice(1) + color_atr*choice(2);

ratio_norm(ratio_norm==0) = NaN;

s(1) = figure(1);
plot(ratio_norm', 'color', clr(1,:), 'linewidth', 2);
ylim(range_overlap);
hold on
plot(nanmean(ratio_norm), 'color', clr(2,:), 'linewidth', 4);
% plot([1 size(ratio_norm,2)], [1 1], ':r', 'linewidth', 1.5);
set(gca, 'visible', 'off')
set(gcf, 'color', 'w')

s(2) = figure(2);
imagesc(ratio_norm); colormap gray; 
caxis(range_separated);
set(gca, 'ticklength', [0 0], 'xticklabel', [], 'yticklabel', []);
c = colorbar; colormap(flipud(gray));
set(c, 'yticklabel', []);

end