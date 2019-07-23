function [h] = plot_waterfall(ratio, colornum)

ratio_smd = reshape(smooth(ratio),[],size(ratio,2));
cmp_wtf = hsv(colornum);

h = waterfall(ratio_smd');
set(h, 'FaceColor', 'flat');
set(h, 'FaceAlpha', 0.5);
set(h, 'EdgeColor', 'k');
set(h, 'FaceVertexCData', cmp_wtf(1:size(ratio, 2), :));

set(gcf, 'Color', [1 1 1]);

set(gca, 'Color', get(gcf, 'Color'));
set(gca, 'GridLineStyle', 'none');

end