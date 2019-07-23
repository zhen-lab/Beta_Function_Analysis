% cmp = lbmap(size(ratio_rev, 2), 'redblue');
cmp = hsv(size(ratio_rev, 2));
fps = 17.22;

figure(1); 
for i = 1 : size(ratio_rev, 2)
subplot(size(ratio_rev, 2), 1 , i); ratio_norm = normalize_signal(ratio_rev(:, size(ratio_rev, 2)-i+1));
plot(ratio_norm, 'color', cmp(i, :), 'linewidth', 2.5); hold on;
curv_norm = normalize_signal(curv_rev(:, size(ratio_rev, 2)-i+1));
plot(curv_norm, 'color', 'k', 'linewidth', 1); box off; ylim([0 1.2]); xlim([1 size(ratio_rev, 1)]);
if i <size(ratio_rev, 2)
set(gca, 'xcolor', [1 1 1]);
end; set(gca, 'ytick', [0 1]); set(gca, 'xtick', 0:floor(size(ratio_rev, 1)/fps)*5:size(ratio_rev, 1));
end; set(gcf, 'color', [1 1 1]); hold off;
set(gca, 'xtick', 0:fps*2:size(ratio_rev, 1));
tick = get(gca, 'xtick');
set(gca, 'xticklabel', tick/fps);

figure(2); 
for i = 1 : size(ratio_rev, 2)-1
subplot(size(ratio_rev, 2), 1 , i); ratio_norm = normalize_signal(ratio_rev(:, size(ratio_rev, 2)-i+1));
plot(ratio_norm, 'color', cmp(i, :), 'linewidth', 2.5); hold on;
curv_norm = normalize_signal(curv_rev(:, size(ratio_rev, 2)-i));
plot(curv_norm, 'color', 'k', 'linewidth', 1); box off; ylim([0 1.2]); xlim([1 size(ratio_rev, 1)]);
if i <size(ratio_rev, 2)-1
set(gca, 'xcolor', [1 1 1]);
end; set(gca, 'ytick', [0 1]); set(gca, 'xtick', 0:size(ratio_rev, 1)/fps:size(ratio_rev, 1));
end; set(gcf, 'color', [1 1 1]); hold off;
set(gca, 'xtick', 0:fps*2:size(ratio_rev, 1));
tick = get(gca, 'xtick');
set(gca, 'xticklabel', tick/fps);

figure(3); 
for i = 1 : size(ratio_rev, 2)-1
subplot(size(ratio_rev, 2), 1 , i+1); ratio_norm = normalize_signal(ratio_rev(:, size(ratio_rev, 2)-i));
plot(ratio_norm, 'color', cmp(i+1, :), 'linewidth', 2.5); hold on;
curv_norm = normalize_signal(curv_rev(:, size(ratio_rev, 2)-i+1));
plot(curv_norm, 'color', 'k', 'linewidth', 1); box off; ylim([0 1.2]); xlim([1 size(ratio_rev, 1)]);
if i <size(ratio_rev, 2)-1
set(gca, 'xcolor', [1 1 1]);
end; set(gca, 'ytick', [0 1]); set(gca, 'xtick', 0:size(ratio_rev, 1)/fps:size(ratio_rev, 1));
end; set(gcf, 'color', [1 1 1]); hold off;
set(gca, 'xtick', 0:fps*2:size(ratio_rev, 1));
tick = get(gca, 'xtick');
set(gca, 'xticklabel', tick/fps);