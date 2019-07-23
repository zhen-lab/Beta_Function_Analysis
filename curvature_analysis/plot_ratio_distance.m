ratio_smd = reshape(smooth(ratio), [], size(ratio,2));
ratio_smd = reshape(smooth(ratio_smd), [], size(ratio,2));
ratio_max = max(ratio_smd);
ratio_max = repmat(ratio_max, size(ratio, 1), 1);
ratio_min = min(ratio);
ratio_min = repmat(ratio_min, size(ratio, 1), 1);
ratio_del = ratio-ratio_min;
ratio_rng = ratio_max-ratio_min;
ratio_norm = ratio_del./ratio_rng;

distSquare_smd = reshape(smooth(distSquare), [], size(distSquare,2));
distSquare = distSquare_smd';
distSquare_max = max(distSquare);
distSquare_max = repmat(distSquare_max, size(distSquare,1), 1);
distSquare_min = repmat(min(distSquare), size(distSquare,1), 1);
distSquare_del = distSquare-distSquare_min;
distSquare_rng = distSquare_max-distSquare_min;
distSquare_norm = distSquare_del./distSquare_rng;

figure(1);

set(gcf, 'Color', [1 1 1]);
cmp=hsv(size(ratio, 2));

for i=size(ratio, 2)-1:-1:1
    
    axes=subplot(size(ratio, 2) - 1, 1, i);
    plot(ratio_norm(:, i), 'Color', cmp(i, :), 'LineWidth', 1); hold on;
    plot(distSquare_norm(:, i), 'Color',  [1 1 1]*0.2,  'LineWidth', 1.5); hold on;
    plot(ratio_norm(:, i+1), 'Color', cmp(i+1, :), 'LineWidth', 1);
    leg=legend(['DA ' num2str(i+1)], 'Distance', ['DA ' num2str(i+2)], 'Location', 'EastOutside'); 
    set(leg, 'Box', 'off'); 
    set(gca, 'XTickLabel', [], 'XColor', get(gcf, 'Color')); ylim([0 1.5]);
    
    if i==size(ratio, 2)-1 
        xlabel('Frame Number'); 
        leg=legend(['DA ' num2str(i+1)], 'Distance', ['DA ' num2str(i+2) '/' num2str(i+3)], 'Location', 'EastOutside'); 
        set(leg, 'Box', 'off'); 
        set(gca, 'XTickLabel', 0:floor(length(ratio_norm)/10):length(ratio_norm), 'XColor', [1 1 1]*0.2);
    end
    
    set(gca, 'Color', get(gcf, 'Color')); 
    box off; 
    set(axes, 'TickLength', [0 0]);

end

title('\color{red}Activity \color{black}and \color{red}Pairwise Distance \color{black}of DA Neurons in \itunc-13 \rmL1 larva');
