function [max_loc] = plot_xcorr_stitched( ratio_rev, curv_rev )

ratio_rev_lin = reshape(ratio_rev, [], 1);
curv_rev_lin = reshape(curv_rev, [], 1);

[a,b] = xcorr(ratio_rev_lin, curv_rev_lin, 'biased');

figure(2);
subplot(3,1,1); plot(ratio_rev_lin, 'Color', [0.9,0.1,0], 'LineWidth', 2); axis tight; ylabel('Activity');
subplot(3,1,2); plot(curv_rev_lin, 'Color', [0,0.5,0.9], 'LineWidth', 2); axis tight; ylabel('Curvature');
subplot(3,1,3); plot(b, a, 'Color', [0.5,0.5,0.5], 'LineWidth', 2); axis tight; ylabel('Cross-correlation');

max_loc = b(a == max(a));

end