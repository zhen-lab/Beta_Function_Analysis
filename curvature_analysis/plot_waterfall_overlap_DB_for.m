function [] = plot_waterfall_overlap_DB_for(ratio, curv, ratio_for, curv_for, colornum)

if size(ratio)==size(curv)
    
    numframe = size(ratio, 1);
    fps = 10;
    
    figure(1); 

    set(gcf, 'Color', [1 1 1]);
    
    az = -5;
    el = 78;

    ratio_smd = reshape(smooth(ratio), [], size(ratio, 2));
    curv_smd = reshape(smooth(curv), [], size(curv, 2));
    ratio_for = reshape(smooth(ratio_for), [], size(ratio_for, 2));
    curv_for = reshape(smooth(curv_for), [], size(curv_for, 2));
    cmp_wtf = hsv(colornum);

    h = waterfall(ratio_smd(:, 1:end-1)');
    hold on;
    s = waterfall(curv_smd(:, 2:end)');
    hr = waterfall(ratio_for(:, 1:end-1)');
    sr = waterfall(curv_for(:, 2:end)');

    set(h, 'FaceColor', 'flat');
    set(h, 'FaceAlpha', 0.2);
    set(h, 'EdgeColor', 'k');
    set(h, 'FaceVertexCData', cmp_wtf(1:colornum-1, :));
    set(hr, 'FaceColor', 'flat');
    set(hr, 'FaceAlpha', 0.5);
    set(hr, 'EdgeColor', 'none');
    set(hr, 'FaceVertexCData', cmp_wtf(1:colornum-1, :));
    
    set(s, 'FaceColor', 'flat');
    set(s, 'FaceAlpha', 0.3);
    set(s, 'EdgeColor', 'k');
    set(s, 'LineStyle', '--');
    set(s, 'FaceVertexCData', cmp_wtf(1:colornum-1, :));
    set(sr, 'FaceColor', 'flat');
    set(sr, 'FaceAlpha', 0.6);
    set(sr, 'EdgeColor', 'none');
    set(sr, 'FaceVertexCData', cmp_wtf(1:colornum-1, :));

    set(gca, 'Color', get(gcf, 'Color'));
    set(gca, 'GridLineStyle', 'none');

    title('\color{red}Activity \color{black}and \color{red}Curvature \color{black}of B motoneurons in locomotion');
    xlabel('Time/s'); ylabel('Neuron number'); zlabel('Signal');
    
    set(gca, 'XTick', 0:floor(numframe/10):numframe);
    xTick = get(gca, 'XTick');
    set(gca, 'XTickLabel', floor(xTick/fps));
    set(gca, 'YTick', 1:size(ratio,2));
    set(gca, 'YTickLabel', {'2/3'; '3/4'; '4/5';'5/6';'6/7'; '8/7'});
    view(az,el);
    axis tight;

    if size(ratio, 2)==2
        
        ylim([0.5 2.5]);
    
    end

else
    
    disp('Sizes do not match.');
    
end

end