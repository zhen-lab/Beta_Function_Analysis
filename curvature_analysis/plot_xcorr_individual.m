function [ coeff, lags ] = plot_xcorr_individual( ratio_rev_3_9_norm, curv_rev_2_7_norm)

dim_1 = size(ratio_rev_3_9_norm, 1);
dim_2 = size(ratio_rev_3_9_norm, 2);

reptime = 10;

coeff = zeros(2 * dim_1 * reptime - 1, dim_2);
lags = zeros(dim_2, 2 * dim_1 * reptime - 1);

cmp = hsv(dim_2);

for iii = 1:size(ratio_rev_3_9_norm, 2)
   
    ratio_rev_3_9_norm_rep = repmat(ratio_rev_3_9_norm(:, iii), reptime, 1);
    curv_rev_2_7_norm_rep = repmat(curv_rev_2_7_norm(:, iii), reptime, 1);
    
    [coeff(:, iii), lags(iii, :)] = xcorr(ratio_rev_3_9_norm_rep, curv_rev_2_7_norm_rep, 'biased');
    plot(lags(iii, :), coeff(:, iii), 'Color', cmp(iii, :), 'LineWidth', 2);
    hold on;
    
end

end