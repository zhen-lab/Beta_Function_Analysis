pos_neg_activity;

figure;

dorsal_neg_body = dorsal_neg(35:end, :);
dorsal_pos_body = dorsal_pos(35:end, :);
ventral_neg_body = ventral_neg(35:end, :);
ventral_pos_body = ventral_pos(35:end, :);

vec_dorsal_neg_body = reshape(dorsal_neg_body(dorsal_neg_body>0), [], 1);
vec_dorsal_pos_body = reshape(dorsal_pos_body(dorsal_pos_body>0), [], 1);
vec_ventral_neg_body = reshape(ventral_neg_body(ventral_neg_body>0), [], 1);
vec_ventral_pos_body = reshape(ventral_pos_body(ventral_pos_body>0), [], 1);

maxix = max([max(vec_dorsal_neg_body) max(vec_dorsal_pos_body) max(vec_ventral_neg_body) max(vec_ventral_pos_body)]);
maxiy = 0.15;
steps = 70;
xbin = 0:ceil(maxix/steps):ceil(maxix/steps)*steps;

subplot(2,2,1); 
[f1,x1] = hist(vec_dorsal_neg_body, xbin);
bar(x1, f1/sum(f1));
h = findobj(gca, 'Type', 'patch');
set(h, 'facecolor', [1 0.2 0], 'facealpha', 0.4, 'edgecolor', 'w'); 
xlim([-ceil(maxix/steps), maxix]); ylim([0 maxiy]);title('Dorsal Contraction', 'fontsize', 20);

subplot(2,2,3);
[f1,x1] = hist(vec_dorsal_pos_body, xbin);
bar(x1, f1/sum(f1));
h = findobj(gca, 'Type', 'patch');
set(h, 'facecolor', [1 0.2 0], 'facealpha', 0.4, 'edgecolor', 'w'); title('Dorsal Extension', 'fontsize', 20);
xlim([-ceil(maxix/steps), maxix]); ylim([0 maxiy]);

subplot(2,2,2); 
[f1,x1] = hist(vec_ventral_neg_body, xbin);
bar(x1, f1/sum(f1));
h = findobj(gca, 'Type', 'patch');
set(h, 'facecolor', [0 0.2 1], 'facealpha', 0.4, 'edgecolor', 'w'); 
xlim([-ceil(maxix/steps), maxix]); ylim([0 maxiy]);title('Ventral Extension', 'fontsize', 20);

subplot(2,2,4); 
[f1,x1] = hist(vec_ventral_pos_body, xbin);
bar(x1, f1/sum(f1));
h = findobj(gca, 'Type', 'patch');
set(h, 'facecolor', [0 0.2 1], 'facealpha', 0.4, 'edgecolor', 'w'); title('Ventral Contraction', 'fontsize', 20);
xlim([-ceil(maxix/steps), maxix]); ylim([0 maxiy]);