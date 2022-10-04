function [ coeff, lags ] = plot_xcorr_individual( activities, curvdata, reptime)

dim_1 = size(activities, 1);
dim_2 = size(activities, 2);

% reptime = 10;

coeff = zeros(2 * dim_1 * reptime - 1, dim_2);
lags = zeros(dim_2, 2 * dim_1 * reptime - 1);

cmp = hsv(dim_2);

for i = 1:size(activities, 2)
   
    activities_norm_rep = repmat(activities(:, i), reptime, 1);
    curvdata_norm_rep = repmat(curvdata(:, i), reptime, 1);
    
    [coeff(:, i), lags(i, :)] = xcorr(activities_norm_rep, curvdata_norm_rep, 'biased');
    plot(lags(i, :), coeff(:, i), 'Color', cmp(i, :), 'LineWidth', 2);
    hold on;
    
end

end