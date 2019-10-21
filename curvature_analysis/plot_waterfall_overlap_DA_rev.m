function [] = plot_waterfall_overlap_DA_rev(ratio, curv, ratio_rev, curv_rev, colornum)

if size(ratio)==size(curv)
    
    numframe = size(ratio, 1);
    numneuron = size(ratio, 2);
    fps = 10;
    
    figure(1); 

    set(gcf, 'Color', [1 1 1]);
    
    az = -5;
    el = 78;

    ratio_smd = reshape(smooth(ratio), [], size(ratio, 2));
    curv_smd = reshape(smooth(curv), [], size(curv, 2));
    ratio_rev = reshape(smooth(ratio_rev), [], size(ratio_rev, 2));
    curv_rev = reshape(smooth(curv_rev), [], size(curv_rev, 2));
    cmp_wtf = hsv(colornum);
% 
%     h = waterfall(ratio_smd(:, 2:end)');
%     hold on;
%     s = waterfall(curv_smd(:, 1:end-1)');
%     hr = waterfall(ratio_rev(:, 2:end)');
%     sr = waterfall(curv_rev(:, 1:end-1)');
    
    h = waterfall(ratio_smd');
    hold on;
    s = waterfall(curv_smd');
    hr = waterfall(ratio_rev');
    sr = waterfall(curv_rev');

    set(h, 'FaceColor', 'flat');
    set(h, 'FaceAlpha', 0.2);
    set(h, 'EdgeColor', 'k');
    set(h, 'FaceVertexCData', cmp_wtf(1:numneuron, :));
    set(hr, 'FaceColor', 'flat');
    set(hr, 'FaceAlpha', 0.5);
    set(hr, 'EdgeColor', 'none');
    set(hr, 'FaceVertexCData', cmp_wtf(1:numneuron, :));
    
    set(s, 'FaceColor', 'flat');
    set(s, 'FaceAlpha', 0.3);
    set(s, 'EdgeColor', 'k');
    set(s, 'LineStyle', '--');
    set(s, 'FaceVertexCData', cmp_wtf(1:numneuron, :));
    set(sr, 'FaceColor', 'flat');
    set(sr, 'FaceAlpha', 0.6);
    set(sr, 'EdgeColor', 'none');
    set(sr, 'FaceVertexCData', cmp_wtf(1:numneuron, :));

    set(gca, 'Color', get(gcf, 'Color'));
    set(gca, 'GridLineStyle', 'none');

    title('\color{red}Activity \color{black}and \color{red}Curvature \color{black}of A motoneurons in locomotion');
    xlabel('Time/s'); ylabel('Neuron number'); zlabel('Signal');
    
    set(gca, 'XTick', 0:floor(numframe/10):numframe);
    xTick = get(gca, 'XTick');
    set(gca, 'XTickLabel', floor(xTick/fps));
    set(gca, 'YTick', 1:size(ratio,2));
    set(gca, 'YTickLabel', {'3/2'; '4/3'; '5/4';'6/5';'7/6'; '8&9/7'});
    view(az,el);
    axis tight;

    if size(ratio, 2)==2
        
        ylim([0.5 2.5]);
    
    end

else
    
    disp('Sizes do not match.');
    
end

end