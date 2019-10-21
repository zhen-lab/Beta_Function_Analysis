function [ ratio_rev_3_9_norm, ratio_rev_lin_3_9_norm, curv_rev_2_7_norm, curv_rev_lin_2_7_norm ] = plot_waterfall_subplot_norm( ratio_rev, curv_rev, colornum )

ratio_rev_3_9 = ratio_rev(:, 2:7);
ratio_rev_3_9_max = max(ratio_rev_3_9);
ratio_rev_3_9_max = repmat(ratio_rev_3_9_max, size(ratio_rev_3_9, 1), 1);
ratio_rev_3_9_min = min(ratio_rev_3_9);
ratio_rev_3_9_min = repmat(ratio_rev_3_9_min, size(ratio_rev_3_9, 1), 1);
ratio_rev_3_9_del = ratio_rev_3_9-ratio_rev_3_9_min;
ratio_rev_3_9_rng = ratio_rev_3_9_max-ratio_rev_3_9_min;
ratio_rev_3_9_norm = ratio_rev_3_9_del./ratio_rev_3_9_rng;
ratio_rev_lin_3_9_norm = reshape(ratio_rev_3_9_norm, [], 1);

curv_rev_2_7 = curv_rev(:, 1:6);
curv_rev_2_7_max = max(curv_rev_2_7);
curv_rev_2_7_max = repmat(curv_rev_2_7_max, size(curv_rev_2_7,1), 1);
curv_rev_2_7_min = repmat(min(curv_rev_2_7), size(curv_rev_2_7,1), 1);
curv_rev_2_7_del = curv_rev_2_7-curv_rev_2_7_min;
curv_rev_2_7_rng = curv_rev_2_7_max-curv_rev_2_7_min;
curv_rev_2_7_norm = curv_rev_2_7_del./curv_rev_2_7_rng;
curv_rev_lin_2_7_norm = reshape(curv_rev_2_7_norm, [], 1);

plot_waterfall_subplot(ratio_rev_3_9_norm, curv_rev_2_7_norm, colornum);

end