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

% fps = 10;
% framenum = size(ratio, 1);
cmp = hsv(size(ratio, 2));
size_shifted = size(ratio, 2)-1;

ratio_smd = reshape(smooth(ratio), [], size(ratio, 2));
curv_smd = reshape(smooth(curv), [], size(curv, 2));

range_ru = max(max(ratio_smd(:, 1:end)));
range_rd = min(min(ratio_smd(:, 1:end)));
range_cu = max(max(curv_smd(:, 2:end)));
range_cd = min(min(curv_smd(:, 2:end)));

for i = size_shifted:-1:1
    
    % Plot ratio
    subplot(size_shifted, 2, 2*i-1);
    ratio_smd = smooth(ratio(:, size(ratio, 2)-i));
    plot(ratio_smd, 'Color', cmp(size(ratio, 2)-i, :), 'LineWidth', 1.5);
    hold on;
    line([boundary(1) boundary(1)], [range_rd range_ru], 'LineStyle', '--', 'Color', [0.2 0.2 0.2]);
    line([boundary(2) boundary(2)], [range_rd range_ru], 'LineStyle', '--', 'Color', [0.2 0.2 0.2]);    
    leg = legend(['DB ' num2str(size(ratio, 2)-i+1)], 'Location', 'NorthEast'); 
    set(leg, 'Box', 'off');
    set(gca, 'Color', get(gcf, 'Color')); 
    box off;
    set(gca, 'XColor', get(gcf, 'Color'));
%     set(gca, 'XTick', 0:floor(framenum/fps):framenum);
%     x_tick = get(gca, 'XTick');
%     set(gca, 'XTickLabel', floor(x_tick/fps));
        
%     if i == size_shifted
%         
% %         xlabel('Time/s');
% %         set(gca, 'XTick', 0:floor(framenum/fps):framenum);
% %         x_tick = get(gca, 'XTick');
% %         set(gca, 'XTickLabel', floor(x_tick/fps));
%         leg=legend(['DB ' num2str(size(ratio, 2)-i+1)], 'Location', 'NorthEastOutside'); 
%         set(leg, 'Box', 'off'); 
%         set(gca, 'XTickLabel', 0:floor(length(ratio)/10):length(ratio), 'XColor', [1 1 1]*0.2);
    
%     end
    
    hold on; 
    ylim([range_rd range_ru]);
%     ylabel('GCaMP6/RFP');
%     xlabel('Time/s');

    % Plot curvature    
    subplot(size_shifted, 2, 2*i);
    curv_smd = smooth(curv(:, size(ratio, 2)-i+1));
    plot(curv_smd, 'Color',  cmp(size(ratio, 2)-i, :),  'LineWidth', 1.5);
    line([boundary(1) boundary(1)], [range_cd range_cu], 'LineStyle', '--', 'Color', [0.2 0.2 0.2]);
    line([boundary(2) boundary(2)], [range_cd range_cu], 'LineStyle', '--', 'Color', [0.2 0.2 0.2]);    
    leg=legend(['DB ' num2str(size(ratio, 2)-i+2)], 'Location', 'NorthEast'); 
    set(leg, 'Box', 'off'); 
    set(gca, 'XTickLabel', [], 'XColor', get(gcf, 'Color'));
    box off;
%     set(gca, 'XTick', 0:floor(framenum/fps):framenum);
%     x_tick = get(gca, 'XTick');
%     set(gca, 'XTickLabel', floor(x_tick/fps));
    set(gca, 'XColor', get(gcf, 'Color'));
    
%     if i==size_shifted
%         
%         xlabel('Time/s');
% %         set(gca, 'XTick', 0:floor(framenum/fps):framenum);
%         x_tick = get(gca, 'XTick');
%         set(gca, 'XTickLabel', floor(x_tick/fps));
%         leg=legend(['DB ' num2str(size(ratio, 2)-i+2)], 'Location', 'NorthEastOutside'); 
%         set(leg, 'Box', 'off'); 
%         set(gca, 'XTickLabel', 0:floor(length(ratio)/10):length(ratio), 'XColor', [1 1 1]*0.2);
%     
%     end
    
    hold on; 
    ylim([range_cd range_cu]);
%     ylabel('Curvature');
%     xlabel('Time/s');
    
    set(gca, 'Color', get(gcf, 'Color')); 
    box off;

end

% legendInfo{num_data} = ['A motoneuron # ' num2str(num_data+1) '/' num2str(num_data+2)];
% hold off;

% title('B motoneuron activity in L1 larva');

end