bodywall = 35;
dorsal_vec_ab = reshape(dorsal_smd_ab(bodywall:end, :), [] ,1);
ventral_vec_ab = reshape(ventral_smd_ab(bodywall:end, :), [], 1);
dorsal_vec_ctrl = reshape(dorsal_smd_ctrl(bodywall:end, :), [] ,1);
ventral_vec_ctrl = reshape(ventral_smd_ctrl(bodywall:end, :), [], 1);

% xbin = 0:100:20000;
% x_lim = 20000;
% y_lim = 15000;
xbin = 0:200:8000;
x_lim = 8000;
y_lim = 0.5;
% xbin = 0:0.01:2;
% x_lim = 2;
% y_lim = 15000;

[f1,x1] = hist(dorsal_vec_ab, xbin);
[f2,x2] = hist(dorsal_vec_ctrl, xbin);

figure;
bar(x1, f1/sum(f1));
hold on;
bar(x2, f2/sum(f2));
h = findobj(gca, 'Type', 'patch');
set(h(1), 'facecolor', [1 0.2 0], 'facealpha', 0.4, 'edgecolor', 'w'); set(h(2), 'facecolor', [1 0.2 0], 'facealpha', 1, 'edgecolor', 'w');
set(gcf, 'color', [1 1 1 ]); set(gca, 'layer', 'top');
box off;
xlabel('Muscle activity'); ylabel('Density'); title('Histogram of muscle activity on the Dorsal side');
% xlim([-0.5 x_lim]);
ylim([0 y_lim]);

[f1,x1] = hist(ventral_vec_ab, xbin);
[f2,x2] = hist(ventral_vec_ctrl, xbin);

figure;
bar(x1, f1/sum(f1));
hold on;
bar(x2, f2/sum(f2));
h = findobj(gca, 'Type', 'patch');
set(h(1), 'facecolor', [0 0.2 1], 'facealpha', 0.4, 'edgecolor', 'w'); set(h(2), 'facecolor', [0 0.2 1], 'facealpha', 1, 'edgecolor', 'w');
set(gcf, 'color', [1 1 1 ]); set(gca, 'layer', 'top');
box off;
xlabel('Muscle activity'); ylabel('Density'); title('Histogram of muscle activity on the Ventral side');
% xlim([-0.5 x_lim]);
ylim([0 y_lim]);

means = [mean(mean(dorsal_smd_ab)); mean(mean(dorsal_smd_ctrl)); mean(mean(ventral_smd_ab)); mean(mean(ventral_smd_ctrl))]
