function [] = plot_heatmap_perievent(ratio)

opengl software;

figure(1);

subplot(2,1,1);
imagesc(ratio'); colormap(hot);
title('Heatmap of DB activities', 'fontsize', 15);
set(gca, 'YDir', 'normal');set(gca, 'ytick', 1:size(ratio, 2)); % To display the head-tail direction from bottom to top
set(gca, 'yticklabel', {'2' '3' '4' '5' '6' '7' '8/9'}); fps = 17.22; numframes = size(ratio, 1);
ylabel('DA neuron number');  set(gca, 'xtick', [], 'ticklength', [0.001 0.001], 'linewidth', 2.5); box on;
xlim([1 size(ratio, 1)+1]);

subplot(2,1,2); 
%fill([1 1 958 958], [0 10 10 0], [0.8 0.8 0.8], 'edgecolor', 'none');
hold on;
cmp = lbmap(size(ratio, 2), 'redblue');
for i=1:size(ratio, 2)
plot(ratio(:, i), 'color', cmp(i, :), 'linewidth', 2);
end;
title('Peri-event plot of DB activities', 'fontsize', 15);
xlim([1 size(ratio, 1)+1]); ylim([0 1.1*max(max(ratio))]);
set(gca, 'xtick', [], 'ticklength', [0.001 0.001], 'linewidth', 2.5); 
set(gca,'XTICK',1:10*fps:numframes);
x_tick=get(gca,'XTICK');
set(gca,'XTICKLABEL',(x_tick-1)/fps, 'ticklength', [0.001 0.001], 'linewidth', 2.5, 'layer', 'top');xlabel('Time/s'); ylabel('Signal');
set(gcf, 'color', [1 1 1 ]);

hold off;

end
