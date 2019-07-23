function [] = plot_waterfall_subplot(ratio_rev, curv_rev, colornum)

figure(1); 

az = -10; % Set up viewing angle
el = 60;
fps = 17.22;
framenum = size(ratio_rev, 1);

subplot(1,2,1); % Plot ratio
plot_waterfall(ratio_rev, colornum);

title('\color{red}Activity \color{black}of A motoneurons');
xlabel('Time/s'); ylabel('Neuron number'); zlabel('GCaMP6/wCherry');
set(gca, 'XTick', 0:floor(framenum/fps):framenum);
x_tick = get(gca, 'XTick');
set(gca, 'XTickLabel', floor(x_tick/fps));
set(gca, 'YTick', 1:size(ratio_rev,2));
set(gca, 'YTickLabel', {2;3;4;5;6;7;'8/9'});
view(az,el);
axis tight;

subplot(1,2,2); % Plot curvature
plot_waterfall(curv_rev, colornum); 

title('\color{red}Curvature \color{black}of A motoneurons');
xlabel('Time/s'); ylabel('Neuron number'); zlabel('Curvature');
set(gca, 'XTick', 0:floor(framenum/fps):framenum);
x_tick = get(gca, 'XTick');
set(gca, 'XTickLabel', floor(x_tick/fps));
set(gca, 'YTick', 1:size(curv_rev,2));
set(gca, 'YTickLabel', {2;3;4;5;6;7;'8/9'});
view(az,el);
axis tight;

end