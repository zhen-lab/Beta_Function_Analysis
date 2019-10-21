function [] = plot_overlap_together(ratio, boundary, colornum)

figure(1);

framenum = size(ratio, 1);
numneurons = size(ratio, 2);
cmp = hsv(colornum);
fps = 10;

for i = 1:numneurons-1
    
    ratio_smd = smooth(ratio(:, i));
    plot(ratio_smd, 'Color', cmp(i, :), 'LineWidth', 2);
    hold on;
    plot(ratio_smd(boundary), 'LineWidth', 5, 'Color', cmp(i, :));
    
end

xlabel('Time/s');
ylabel('GCaMP6/RFP');
title('L1 DB motorneuron activity change in locomotion');

set(gca, 'XTick', 0:floor(framenum/fps):framenum);
x_tick = get(gca, 'XTick');
set(gca, 'XTickLabel', floor(x_tick/fps));

end