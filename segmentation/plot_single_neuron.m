function [] = plot_single_neuron(gcamp, rfp, nm)

figure; 

subplot(211); hold on
plot(rfp(:, nm), 'color', [1 0.3 0], 'linewidth', 3);
plot(gcamp(:, nm), 'color', [0 1 0], 'linewidth', 3);
ylimit = get(gca, 'ylim');
for i = 0:4
    plot([100+(i-1)*100 100+(i-1)*100], ylimit, '--', 'color', 'k', 'linewidth', 3)
end
box off;
set(gca, 'ticklength', [0 0], 'xticklabels', []);

subplot(212); hold on
plot(gcamp(:, nm)./rfp(:, nm), 'k', 'linewidth', 3);
ylimit = get(gca, 'ylim');
for i = 0:4
    plot([100+(i-1)*100 100+(i-1)*100], ylimit, '--', 'color', 'k', 'linewidth', 3)
end
box off;
set(gca, 'ticklength', [0 0], 'xticklabels', []);

end
