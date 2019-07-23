function [] = plot_xcorr_separated( ratio_rev, curv_rev )

a = zeros(size(ratio_rev, 2), 2*size(ratio_rev, 1)-1);
b = zeros(2*size(ratio_rev, 1)-1, size(ratio_rev, 2));

for n_num = 1:size(ratio_rev, 2)
   
    [a(n_num, :), b(:, n_num)] = xcorr(ratio_rev(:, n_num), curv_rev(:, n_num), 'none');
    
end

figure(3);

plot_waterfall(a');

title('Cross-correlation between activity and curvature of A motoneurons in L1 during backward movement');
xlabel('Frame number'); ylabel('Neuron number'); zlabel('Cross-correlation');

b_lim = b(:, 1); % Coordinates for x axis
mpt = ceil(length(b_lim)/2);
step = ceil(length(b_lim)/20);
set(gca, 'xtick', mpt-ceil(mpt/step)*step:step:mpt+ceil(mpt/step)*step);
set(gca, 'xticklabel', -ceil(mpt/step)*step:step:ceil(mpt/step)*step); % This makes sure point zero is labeled
axis tight;

end