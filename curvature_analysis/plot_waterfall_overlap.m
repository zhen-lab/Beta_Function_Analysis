function [] = plot_waterfall_overlap(ratio, curv, colornum)

if size(ratio)==size(curv)
    
    numframe = size(ratio, 1);
    fps = 17;
    
    figure(1); 

    set(gcf, 'Color', [1 1 1]);
    
    az = -5;
    el = 70;

    ratio_smd = reshape(smooth(ratio), [], size(ratio, 2));
    curv_smd = reshape(smooth(curv), [], size(curv, 2));
    cmp_wtf = hsv(colornum);

    h = waterfall(ratio_smd(:, 2:end)');
    hold on;
    s = waterfall(curv_smd(:, 1:end-1)');

    set(h, 'FaceColor', 'flat');
    set(h, 'FaceAlpha', 0.45);
    set(h, 'EdgeColor', 'k');
    set(h, 'FaceVertexCData', cmp_wtf(1:size(ratio, 2)-1, :));

    set(s, 'FaceColor', 'flat');
    set(s, 'FaceAlpha', 0.45);
    set(s, 'EdgeColor', 'k');
    set(s, 'LineStyle', '--');
    set(s, 'FaceVertexCData', cmp_wtf(1:size(curv, 2)-1, :));

    set(gca, 'Color', get(gcf, 'Color'));
    set(gca, 'GridLineStyle', 'none');

    title('\color{red}Activity \color{black}and \color{red}Curvature \color{black}of A motoneurons in backward movement');
    xlabel('Time/s'); ylabel('Neuron number'); zlabel('Activity and Curvature');
    
    set(gca, 'XTick', 0:floor(numframe/10):numframe);
    xTick = get(gca, 'XTick');
    set(gca, 'XTickLabel', floor(xTick/fps));
    set(gca, 'YTick', 1:size(ratio,2));
    set(gca, 'YTickLabel', {'2/3';'3/4';'4/5';'5/6';'6/7'});
    view(az,el);
    axis tight;

    if size(ratio, 2)==2
        
        ylim([0.5 2.5]);
    
    end

else
    
    disp('Sizes do not match.');
    
end

end