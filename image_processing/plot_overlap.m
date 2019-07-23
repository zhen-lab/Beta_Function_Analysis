function [] = plot_overlap_separate(ratio, curv, boundary)

% numframes = size(ratio, 1);
% numneurons = size(ratio, 2);
% cmp = hsv(numneurons-1);
% legendInfo = cell(1, numneurons-1);
% fps = 10;

% for i=1:numneurons-1
%     
%     subplot(numneurons-1, 2, 2*i-1);
%     plot(ratio(:,i), 'LineWidth', 2, 'Color', cmp(i, :));
%     
%     set(gca,'XTICK',1:floor(numframes/fps):numframes);
%     x_tick=get(gca,'XTICK');
%     set(gca,'XTICKLABEL', floor(x_tick/fps));
%     xlabel('Time/s');
%     ylabel('GFP/RFP');
%     
%     subplot(numneurons-1, 2, 2*i);
%     plot(curv(:, i+1), 'LineWidth', 2, 'Color', cmp(i, :));
%     
%     set(gca,'XTICK',1:floor(numframes/fps):numframes);
%     x_tick=get(gca,'XTICK');
%     set(gca,'XTICKLABEL', floor(x_tick/fps));
%     xlabel('Time/s');
%     ylabel('Curvature');
%     
%     legendInfo{i} = ['B motoneuron # ' num2str(i+1)];
% %     legend(legendInfo, 'Location', 'BestOutside');
% 
% end;

figure(1);

set(gcf, 'Color', [1 1 1]);
cmp=hsv(size(ratio, 2));
range_ru = max(max(ratio(:, 1:end-1)));
range_rd = min(min(ratio(:, 1:end-1)));
range_cu = max(max(curv(:, 2:end)));
range_cd = min(min(curv(:, 2:end)));

for i = size(ratio, 2)-1:-1:1
    
    % Plot ratio
    axes = subplot(size(ratio, 2) - 1, 2, 2*i-1);
    plot(ratio(:, i), 'Color', cmp(i, :), 'LineWidth', 1.5);
    leg=legend(['DB ' num2str(i+1)], 'Location', 'NorthEast'); 
    set(leg, 'Box', 'off');
    set(gca, 'Color', get(gcf, 'Color')); 
    box off;
    set(axes, 'TickLength', [0 0]);
    set(gca, 'XColor', get(gcf, 'Color'));
    
    if i==size(ratio, 2)-1
        
        xlabel('Frame Number');
        leg=legend(['DB ' num2str(i+1)], 'Location', 'NorthEast'); 
        set(leg, 'Box', 'off'); 
        set(gca, 'XTickLabel', 0:floor(length(ratio)/10):length(ratio), 'XColor', [1 1 1]*0.2);
    
    end
    
    hold on; area(ratio(boundary, i), range_rd, 'FaceColor', cmp(i, :), 'EdgeColor', 'none'); alpha(.3);
    ylim([range_rd range_ru]);
    ylabel('GCaMP6/RFP');

    % Plot curvature    
    axes = subplot(size(ratio, 2) - 1, 2, 2*i);
    plot(curv(:, i+1), 'Color',  cmp(i, :),  'LineWidth', 1.5);
    leg=legend(['DB ' num2str(i+2)], 'Location', 'NorthEast'); 
    set(leg, 'Box', 'off'); 
    set(gca, 'XTickLabel', [], 'XColor', get(gcf, 'Color'));
    box off;
    set(axes, 'TickLength', [0 0]);
    set(gca, 'XColor', get(gcf, 'Color'));
    
    if i==size(ratio, 2)-1
        
        xlabel('Frame Number');
        leg=legend(['DB ' num2str(i+2)], 'Location', 'NorthEast'); 
        set(leg, 'Box', 'off'); 
        set(gca, 'XTickLabel', 0:floor(length(ratio)/10):length(ratio), 'XColor', [1 1 1]*0.2);
    
    end
    
    hold on; area(curv(boundary, i+1), range_cd, 'FaceColor', cmp(i, :), 'EdgeColor', 'none'); alpha(.3);
    ylim([range_cd range_cu]);
    ylabel('Curvature');

    set(gca, 'Color', get(gcf, 'Color')); 
    box off; 
    set(axes, 'TickLength', [0 0]);

end

% legendInfo{num_data} = ['A motoneuron # ' num2str(num_data+1) '/' num2str(num_data+2)];
% hold off;

% title('B motoneuron activity in L1 larva');

end